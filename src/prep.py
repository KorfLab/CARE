# prep.py

import argparse
import os
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


############
# argparse #
############

parser = argparse.ArgumentParser(
	description="Batch subset genome using percentages from YAML")
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


########
# main #
########

outdir = args.output
os.makedirs(args.output, exist_ok=True)

gabs = os.path.abspath(args.genome)
gdir = os.path.dirname(args.genome)
gbsn = os.path.splitext(os.path.basename(gabs))[0]

with open(args.yaml) as ymlin:
	spec = yaml.safe_load(ymlin)

pcts = spec.get("genome_pcts", [])

pct_genome_file_path = {}

for p in pcts:
	if not isinstance(p, float) or not (0 < p <= 1):
		print(f"[prep] Invalid percentage in YAML: {p}, must be a float in (0, 1]")
		continue

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


###########
# cleanup #
###########

if args.cleanup:
	print("\n[prep] Cleaning up...")
	toolbox.only_keep_ext(args.output, ext=[".fa", ".minifq.fastq"])
	print("[prep] Files retained after cleanup:")
	for item in sorted(os.listdir(args.output)):
		print("  -", item)

print(f"\n[prep] var-gn preparation complete")
print(f"[prep] All outputs in: {args.output}")
