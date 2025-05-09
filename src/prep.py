# prep.py

import argparse
import os
import shutil
import sys
import yaml

import toolbox
import aligners


ALIGNER_CLASS_MAP = {
	"star": aligners.StarAligner,
	"hisat2": aligners.Hisat2Aligner,
	"bowtie2": aligners.Bowtie2Aligner,
	"bwa": aligners.BwaAligner,
	"minimap2": aligners.Minimap2Aligner
}


#########
# funcs #
#########

def get_read_names(fq_path):
	"""Get read names from a fastq file"""
	names = set()
	with toolbox.smart_open_read(fq_path) as f:
		for i, line in enumerate(f):
			if i % 4 == 0:
				name = line.strip().split()[0][1:]
				names.add(name)
	return names


def extract_reads(fq_path, shared_names, output_path):
	"""Extract reads by read names in a set to fastq"""
	with toolbox.smart_open_read(fq_path) as fin, open(output_path, "w") as fout:
		while True:
			try:
				block = [next(fin) for _ in range(4)]
			except StopIteration:
				break
			name = block[0].strip().split()[0][1:]
			if name in shared_names:
				fout.writelines(block)


def cleanup(dir):
	for item in os.listdir(dir):
		if (
			item.endswith("p.fa") or
			item.endswith("p.fastq") or
			item.endswith("rl.fastq") or
			item.endswith("rc.fastq")
		):
			continue

		path = os.path.join(dir, item)

		if os.path.isdir(path):
			shutil.rmtree(path)
			print(f"[cleanup] Removed directory: {path}")
		else:
			os.remove(path)
			print(f"[cleanup] Removed file: {path}")


############
# argparse #
############

parser = argparse.ArgumentParser(description="Batch subset genome using percentages from YAML")
parser.add_argument("-g", "--genome", required=True,
	help="Genome FASTA file")
parser.add_argument("-y", "--yaml", required=True,
	help="YAML file of various prep specifications")
parser.add_argument("-a", "--aligner", nargs='+',
	help="Aligners to prepare genome folders for")
parser.add_argument("--r1", required=True,
	help="FASTQ file R1")
parser.add_argument("--r2",
	help="FASTQ file R2 (if paired-end)")
parser.add_argument("-n", "--numReads", type=int, default=10000,
	help="Number of reads to keep in fastq (default=10000)")
parser.add_argument("-t", "--thread", type=int, default=4,
	help="Threads for aligners (default=4)")
parser.add_argument("-o", "--output", required=True,
	help="Directory to store all prep outputs, including genome subsets and alignments")
parser.add_argument("-c", "--cleanup", action="store_true",
	help="Remove all aligner directories and intermediary files")
args = parser.parse_args()

for aligner in args.aligner:
	if aligner not in ALIGNER_CLASS_MAP:
		print(f"[prep] ERROR: Unknown aligner '{aligner}'")
		sys.exit(1)


#########
# check #
#########

if not os.path.isfile(args.genome):
	print(f"[prep] ERROR: Genome file not found: {args.genome}")
	sys.exit(1)

if not os.path.isfile(args.yaml):
	print(f"[prep] ERROR: YAML file not found: {args.yaml}")
	sys.exit(1)

if not os.path.isfile(args.r1):
	print(f"[prep] ERROR: FASTQ R1 file not found: {args.r1}")
	sys.exit(1)

if args.r2 and not os.path.isfile(args.r2):
	print(f"[prep] ERROR: FASTQ R2 file not found: {args.r2}")
	sys.exit(1)


##########
# var-gn #
##########

print("[prep] Generating genome variants for var-gn")

outdir = args.output
os.makedirs(args.output, exist_ok=True)

gabs = os.path.abspath(args.genome)
gdir = os.path.dirname(args.genome)
gbsn = os.path.splitext(os.path.basename(gabs))[0]

with open(args.yaml) as ymlin:
	spec = yaml.safe_load(ymlin)

pcts = spec.get("genome_pcts", [])
if not pcts:
	print("[prep] ERROR: No percentages for genome sizes found in YAML, var-gn prep failed")
	sys.exit(1)

pct_genome_file_path = {}

for p in pcts:
	if not isinstance(p, float) or not (0 < p <= 1):
		print(f"[prep] ERROR: Invalid percentage in YAML: {p}, must be a float in (0, 1]")
		sys.exit(1)

	if p == 1:
		dst = os.path.join(outdir, f"{gbsn}-100p.fa")
		toolbox.cp(gabs, dst)
		pct_genome_file_path[p] = dst
		continue

	cmd = [
		"python3", "subset.py", "xfa",
		"-g", args.genome,
		"-p", f"{p:.2f}",
		"-v"
	]
	toolbox.run(cmd)

	src = os.path.join(gdir, f"{gbsn}.shrunk.fa")
	dst = os.path.join(outdir, f"{gbsn}-{int(p*100)}p.fa")
	toolbox.mv(src, dst)
	pct_genome_file_path[p] = dst

min_pct = min(pcts)
min_pct_genome_file = pct_genome_file_path[min_pct]


###################
# index align xfq #
###################

shrunk_fastqs = {}
read_sets = {}

for aligner in args.aligner:
	aligner_dir = os.path.join(args.output, aligner)
	os.makedirs(aligner_dir, exist_ok=True)

	print(f"\n[prep] Processing aligner: {aligner}")

	aln = ALIGNER_CLASS_MAP[aligner](
		genome=min_pct_genome_file,
		outdir=aligner_dir,
		r1=args.r1,
		r2=args.r2,
		threads=args.thread,
		time=False
	)

	aln.index()
	aln.align()
	sam_path = aln.get_sam_path()

	shrunk_path = os.path.join(aligner_dir, f"{aligner}.shrunk.fastq")

	print(f"\n[xfq] Filtering reads for {aligner} referencing SAM")

	cmd = [
		"python3", "subset.py", "xfq",
		"-g", aln.ref,
		"-s", sam_path,
		"-p", "1",
		"--r1", args.r1
	]

	if args.r2:
		cmd += ["--r2", args.r2]

	cmd += ["-v"]
	toolbox.run(cmd)

	r1_base = os.path.basename(args.r1)
	r1_root = os.path.splitext(r1_base)[0]
	r1_shrunk = os.path.join(os.path.dirname(args.r1), f"{r1_root}.shrunk.fastq")
	r1_dst = os.path.join(aligner_dir, f"{r1_root}.shrunk.{aligner}.fastq")
	toolbox.mv(r1_shrunk, r1_dst)

	shrunk_fastqs[aligner] = {}
	shrunk_fastqs[aligner]["r1"] = r1_dst
	read_sets[aligner] = get_read_names(r1_dst)

	if args.r2:
		r2_base = os.path.basename(args.r2)
		r2_root = os.path.splitext(r2_base)[0]
		r2_shrunk = os.path.join(os.path.dirname(args.r2), f"{r2_root}.shrunk.fastq")
		r2_dst = os.path.join(aligner_dir, f"{r2_root}.shrunk.{aligner}.fastq")
		toolbox.mv(r2_shrunk, r2_dst)
		shrunk_fastqs[aligner]["r2"] = r2_dst


#############
# intersect #
#############

print("\n[prep] Intersecting reads across all aligners")

shared_reads = set.intersection(*read_sets.values())
if not shared_reads:
	print("[prep] ERROR: No shared reads found across aligners")
	sys.exit(1)
print(f"[prep] Found {len(shared_reads)} shared reads across all aligners")

smallest_aligner = min(
	shrunk_fastqs,
	key=lambda a: os.path.getsize(shrunk_fastqs[a]["r1"])
)

r1_for_extract = shrunk_fastqs[smallest_aligner]["r1"]
shared_r1 = os.path.join(outdir, "shared_1.fastq")
extract_reads(r1_for_extract, shared_reads, shared_r1)
print(f"[prep] Shared r1 written to {shared_r1}")

if args.r2:
	r2_for_extract = shrunk_fastqs[smallest_aligner]["r2"]
	shared_r2 = os.path.join(outdir, "shared_2.fastq")
	extract_reads(r2_for_extract, shared_reads, shared_r2)
	print(f"[prep] Shared r2 written to {shared_r2}")


##########
# minifq #
##########

cmd = [
	"python3", "minifq.py",
	"--r1", shared_r1,
	"-n", str(args.numReads),
	"-s", "1",
	"-o", outdir,
	"--sort",
	"-v"
]

if args.r2:
	cmd += ["--r2", shared_r2]

toolbox.run(cmd)

shared_r1_minifq = os.path.join(outdir, "shared_1.minifq.fastq")
min_pct_r1 = os.path.join(outdir, f"shared_1.{int(min_pct*100)}p.fastq")
toolbox.mv(shared_r1_minifq, min_pct_r1)

if args.r2:
	shared_r2_minifq = os.path.join(outdir, "shared_2.minifq.fastq")
	min_pct_r2 = os.path.join(outdir, f"shared_2.{int(min_pct*100)}p.fastq")
	toolbox.mv(shared_r2_minifq, min_pct_r2)

print("[prep] var-gn preparation complete")


##########
# var-rl #
##########

print("\n[prep] Generating read-length variants for var-rl")

read_lengths = spec.get("read_lengths", [])
if not read_lengths:
	print("[prep] ERROR: No read lengths found in YAML, var-rl prep failed")
	sys.exit(1)

for rl in read_lengths:
	if not isinstance(rl, int) or not (0 < rl <= 151):
		print(f"[prep] ERROR: Invalid read length in YAML: {rl}, must be an int in (0, 151]")
		sys.exit(1)

	cmd = [
		"python3", "subset.py", "xrl",
		"-k", str(rl),
		"--r1", min_pct_r1
	]

	if args.r2:
		cmd += ["--r2", min_pct_r2]

	cmd.append("-v")
	toolbox.run(cmd)

print("[prep] var-rl preparation complete")


##########
# var-rc #
##########

print("\n[prep] Generating read subsets for var-rc")

read_counts = spec.get("read_counts", [])
if not read_counts:
	print("[prep] ERROR: No read counts found in YAML, var-rc prep failed")
	sys.exit(1)

for rc in read_counts:
	if not isinstance(rc, int) or rc <= 0:
		print(f"[prep] Error: Invalid read count in YAML: {rc}, must be a positive int")
		sys.exit(1)

	print(f"\n[prep] Targeting {rc} reads")

	shared_num_read = toolbox.get_num_reads(shared_r1)

	print(f"[prep] Available reads in shared FASTQ: {shared_num_read}")

	if shared_num_read >= rc:
		print(f"[prep] Enough reads available - sampling {rc} using minifq.py")

		cmd = [
			"python3", "minifq.py",
			"--r1", shared_r1,
			"-n", str(rc),
			"-o", outdir,
			"-s", "1",
			"--sort",
			"-v"
		]
		if args.r2:
			cmd += ["--r2", shared_r2]

		toolbox.run(cmd)

		base_r1 = os.path.splitext(os.path.basename(shared_r1))[0]
		minifq_r1 = os.path.join(outdir, f"{base_r1}.minifq.fastq")
		rc_r1 = os.path.join(outdir, f"{base_r1}.{rc}rc.fastq")
		toolbox.mv(minifq_r1, rc_r1)

		if args.r2:
			base_r2 = os.path.splitext(os.path.basename(shared_r2))[0]
			minifq_r2 = os.path.join(outdir, f"{base_r2}.minifq.fastq")
			rc_r2 = os.path.join(outdir, f"{base_r2}.{rc}rc.fastq")
			toolbox.mv(minifq_r2, rc_r2)

	else:
		print(f"[prep] Not enough reads - extending to {rc} using weaver.py reuse")

		cmd = [
			"python3", "weaver.py", "reuse",
			"--r1", shared_r1,
			"-n", str(rc),
			"-o", outdir,
			"-s", "1",
			"-v"
		]
		if args.r2:
			cmd += ["--r2", shared_r2]

		toolbox.run(cmd)

		base_r1 = os.path.splitext(os.path.basename(shared_r1))[0]
		weaver_r1 = os.path.join(outdir, f"{base_r1}.weaver.fastq")
		rc_r1 = os.path.join(outdir, f"{base_r1}.{rc}rc.fastq")
		toolbox.mv(weaver_r1, rc_r1)

		if args.r2:
			base_r2 = os.path.splitext(os.path.basename(shared_r2))[0]
			weaver_r2 = os.path.join(outdir, f"{base_r2}.weaver.fastq")
			rc_r2 = os.path.join(outdir, f"{base_r2}.{rc}rc.fastq")
			toolbox.mv(weaver_r2, rc_r2)

print("[prep] var-rc preparation complete")


###########
# cleanup #
###########

if args.cleanup:
	print("\n[prep] Cleaning up...")
	cleanup(outdir)
	print("[prep] Files retained after cleanup:")
	for item in sorted(os.listdir(outdir)):
		print("\t-", item)

print(f"\n[prep] CARE prep complete")
print(f"[prep] All outputs in: {outdir}")


###############
# run-care.sh #
###############

print("\n[prep] Generating run-care.sh ")

min_pct_int = int(min_pct * 100)

astr = " ".join(f"\"{aln}\"" for aln in args.aligner)

with open("run-care.sh", "w") as fout:
	fout.write("#!/bin/bash\n")
	fout.write("set -e\n\n")
	fout.write(f'PREP_DIR="{outdir}"\n')
	fout.write("OUT_DIR=\"results\"\n")
	fout.write("THREADS=4\n")
	fout.write(f"ALIGNERS=({astr})\n\n")
	fout.write("mkdir -p \"$OUT_DIR\"\n\n")

	# var-gn
	fout.write("# Experiment 1: var-gn\n")
	fout.write("python3 care.py var-gn \\\n")
	fout.write(f"  --r1 \"$PREP_DIR/shared_1.{min_pct_int}p.fastq\" \\\n")
	if args.r2:
		fout.write(f"  --r2 \"$PREP_DIR/shared_2.{min_pct_int}p.fastq\" \\\n")
	fout.write("  -g \"$PREP_DIR\"/*p.fa \\\n")
	fout.write("  -t $THREADS \\\n")
	fout.write("  -a \"${ALIGNERS[@]}\" \\\n")
	fout.write("  -o \"$OUT_DIR/var-gn\"\n\n")

	# var-rl
	fout.write("# Experiment 2: var-rl\n")
	fout.write("python3 care.py var-rl \\\n")
	fout.write(f"  --r1 \"$PREP_DIR/shared_1.{min_pct_int}p.\"*rl.fastq \\\n")
	if args.r2:
		fout.write(f"  --r2 \"$PREP_DIR/shared_2.{min_pct_int}p.\"*rl.fastq \\\n")
	fout.write(f"  -g \"$PREP_DIR\"/*{min_pct_int}p.fa \\\n")
	fout.write("  -t $THREADS \\\n")
	fout.write("  -a \"${ALIGNERS[@]}\" \\\n")
	fout.write("  -o \"$OUT_DIR/var-rl\"\n\n")

	# var-cpu
	fout.write("# Experiment 3: var-cpu\n")
	fout.write("python3 care.py var-cpu \\\n")
	fout.write(f"  --r1 \"$PREP_DIR/shared_1.{min_pct_int}p.fastq\" \\\n")
	if args.r2:
		fout.write(f"  --r2 \"$PREP_DIR/shared_2.{min_pct_int}p.fastq\" \\\n")
	fout.write(f"  -g \"$PREP_DIR\"/*{min_pct_int}p.fa \\\n")
	fout.write("  -t 1 2 3 4 \\\n")
	fout.write("  -a \"${ALIGNERS[@]}\" \\\n")
	fout.write("  -o \"$OUT_DIR/var-cpu\"\n\n")

	# var-rc
	fout.write("# Experiment 4: var-rc\n")
	fout.write("python3 care.py var-rc \\\n")
	fout.write("  --r1 \"$PREP_DIR\"/shared_1.*rc.fastq \\\n")
	if args.r2:
		fout.write("  --r2 \"$PREP_DIR\"/shared_2.*rc.fastq \\\n")
	fout.write(f"  -g \"$PREP_DIR\"/*{min_pct_int}p.fa \\\n")
	fout.write("  -t $THREADS \\\n")
	fout.write("  -a \"${ALIGNERS[@]}\" \\\n")
	fout.write("  -o \"$OUT_DIR/var-rc\"\n")

print(f"[prep] Generated run-care.sh")
