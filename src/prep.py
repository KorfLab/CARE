# prep.py

import argparse
import os
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


############
# argparse #
############

parser = argparse.ArgumentParser(
	description="Batch subset genome using percentages from YAML")
parser.add_argument("-g", "--genome", required=True, 
	help="Genome FASTA file")
parser.add_argument("-y", "--yaml", required=True,
	help="YAML file of a list of percentages")
parser.add_argument("-a", "--aligner", nargs='+',
	help="Aligners to prepare genome folders for")
parser.add_argument("--r1", required=True, 
	help="FASTQ file R1")
parser.add_argument("--r2",
	help="FASTQ file R2 (if paired-end)")
parser.add_argument("-n", "----numReads", type=int, default=10000,
    help="Number of reads to keep in fastq (default=10000)")
parser.add_argument("-t", "--thread", type=int, default=4,
	help="Threads for aligners (default=4)")
parser.add_argument("-o", "--output", required=True,
	help="Directory to store all prep outputs, including genome subsets and alignments")


args = parser.parse_args()


########
# main #
########

outdir = args.output
os.makedirs(args.output, exist_ok=True)

gabs = os.path.abspath(args.genome)
gdir = os.path.dirname(args.genome)
gbsn = os.path.splitext(os.path.basename(gabs))[0]

with open(args.yaml) as ymlin:
	pcts = yaml.safe_load(ymlin)

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


#######################
# index align and xfq #
#######################

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

	if args.r2:
		r2_base = os.path.basename(args.r2)
		r2_root = os.path.splitext(r2_base)[0]
		r2_shrunk = os.path.join(os.path.dirname(args.r2), f"{r2_root}.shrunk.fastq")
		r2_dst = os.path.join(aligner_dir, f"{r2_root}.shrunk.{aligner}.fastq")
		toolbox.mv(r2_shrunk, r2_dst)


# TODO add shared reads filtering


# Left to refactor
for aligner in args.aligner:

	r1_shrunk = os.path.join(aligner_dir, f"{os.path.splitext(os.path.basename(args.r1))[0]}.shrunk.{aligner}.fastq")
	r1_mini = os.path.join(aligner_dir, f"{os.path.splitext(os.path.basename(args.r1))[0]}.shrunk.{aligner}.minifq.fastq")

	if args.r2:
		r2_shrunk = os.path.join(aligner_dir, f"{os.path.splitext(os.path.basename(args.r2))[0]}.shrunk.{aligner}.fastq")
		r2_mini = os.path.join(aligner_dir, f"{os.path.splitext(os.path.basename(args.r2))[0]}.shrunk.{aligner}.minifq.fastq")

		cmd = f"python3 minifq.py --r1 {r1_shrunk} --r2 {r2_shrunk} -n {args.numReads} -s 1 --sort -v"
	else:
		cmd = f"python3 minifq.py --r1 {r1_shrunk} -n {args.numReads} -s 1 --sort -v"

	print(f"\n[minifq] Sampling {args.numReads} reads for {aligner}")
	run_cmd(cmd)

	if not os.path.exists(r1_mini):
		moved_r1 = os.path.basename(r1_mini)
		run_cmd(f"mv -f {moved_r1} {r1_mini}")

	if args.r2:
		if not os.path.exists(r2_mini):
			moved_r2 = os.path.basename(r2_mini)
			run_cmd(f"mv -f {moved_r2} {r2_mini}")

print(f"\n[INFO] Obtained {int(min_pct * 100)}% reads for aligners: {args.aligner}")
