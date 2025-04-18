# prep.py

import argparse
import os
import yaml

import toolbox


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


args = parser.parse_args()


########
# main #
########

gabs = os.path.abspath(args.genome)
gdir = os.path.dirname(gabs)
gbsn = os.path.splitext(os.path.basename(gabs))[0]

with open(args.yaml) as ymlin:
	pcts = yaml.safe_load(ymlin)

for p in pcts:
	if not isinstance(p, float) or not (0 < p <= 1):
		print(f"Invalid percentage in YAML: {p}, must be a float in (0, 1]")
		continue

	if p == 1:
		dst = os.path.join(gdir, "ref-100p.fa")
		cmd = ["cp", "-f", gabs, dst]
		run_cmd(' '.join(cmd))
		continue

	cmd = [
		"python3", "subset.py", "xfa",
		"-g", args.genome,
		"-p", f"{p:.2f}",
		"-v"
	]
	run_cmd(' '.join(cmd))

	src = os.path.join(gdir, f"{gbsn}.shrunk.fa")
	dst = os.path.join(gdir, f"ref-{int(p * 100)}p.fa")
	cmd = ["mv", "-f", src, dst]
	run_cmd(' '.join(cmd))

min_pct = min(pcts)
min_pct_str = f"{int(min_pct * 100)}p"
min_pct_genome_name = f"ref-{min_pct_str}"
min_pct_genome_file = f"{min_pct_genome_name}.fa"


########
# move #
########

if args.aligner:
	print(f"\n[INFO] Organizing genome subsets for aligners: {args.aligner}")

	refs = [f for f in os.listdir(gdir) if f.startswith("ref-") and f.endswith("p.fa")]

	for aligner in args.aligner:
		aligner_dir = os.path.join(gdir, f"genomes-{aligner}")
		os.makedirs(aligner_dir, exist_ok=True)

		for ref in refs:
			src = os.path.join(gdir, ref)
			dst = os.path.join(aligner_dir, ref)
			cmd = ["cp", "-f", src, dst]
			run_cmd(' '.join(cmd))


#######################
# index align and xfq #
#######################

for aligner in args.aligner:
	aligner_dir = os.path.join(gdir, f"genomes-{aligner}")
	genome_fa = os.path.join(aligner_dir, min_pct_genome_file)

	# indexing
	if aligner == "star":
		tmp_dir = os.path.join(aligner_dir, "tmp")
		os.makedirs(tmp_dir, exist_ok=True)

		msg = run_cmd_capture(f"STAR --runMode genomeGenerate --runThreadN {args.thread} --genomeDir {tmp_dir} --genomeFastaFiles {genome_fa} 2>&1")
		recommended = None
		for line in msg.splitlines():
			if "genomeSAindexNbases" in line and "recommended" in line:
				recommended = line.split()[-1]
				break

		if recommended is not None:
			print(f"[STAR] Using recommended genomeSAindexNbases = {recommended}")
			run_cmd(f"rm -rf {tmp_dir}")
			run_cmd(f"STAR --runMode genomeGenerate --runThreadN {args.thread} --genomeDir {aligner_dir} --genomeFastaFiles {genome_fa} --genomeSAindexNbases {recommended}")
		else:
			print("[STAR] No recommendation found, using default genomeSAindexNbases")
			run_cmd(f"mv {tmp_dir}/* {aligner_dir}")
			run_cmd(f"rm -rf {tmp_dir}")
	elif aligner == "hisat2":
		run_cmd(f"hisat2-build -p {args.thread} {genome_fa} {genome_fa}")
	elif aligner == "bowtie2":
		run_cmd(f"bowtie2-build --threads {args.thread} {genome_fa} {genome_fa}")
	elif aligner == "bwa":
		run_cmd(f"bwa index {genome_fa}")
	elif aligner == "minimap2":
		run_cmd(f"minimap2 -d {os.path.join(aligner_dir, 'ref.mmi')} {genome_fa}")

	# aligning
	sam_out = os.path.join(aligner_dir, f"subset-{min_pct_genome_name}.sam")
	if aligner == "star":
		if args.r2:
			cmd = f"STAR --runMode alignReads --runThreadN {args.thread} --genomeDir {aligner_dir} --readFilesIn {args.r1} {args.r2} --outFileNamePrefix {aligner_dir}/subset- --outSAMtype SAM"
		else:
			cmd = f"STAR --runMode alignReads --runThreadN {args.thread} --genomeDir {aligner_dir} --readFilesIn {args.r1} --outFileNamePrefix {aligner_dir}/subset- --outSAMtype SAM"
	elif aligner == "hisat2":
		if args.r2:
			cmd = f"hisat2 -p {args.thread} -x {genome_fa} -1 {args.r1} -2 {args.r2} -S {sam_out}"
		else:
			cmd = f"hisat2 -p {args.thread} -x {genome_fa} -U {args.r1} -S {sam_out}"
	elif aligner == "bowtie2":
		if args.r2:
			cmd = f"bowtie2 -p {args.thread} -x {genome_fa} -1 {args.r1} -2 {args.r2} -S {sam_out}"
		else:
			cmd = f"bowtie2 -p {args.thread} -x {genome_fa} -U {args.r1} -S {sam_out}"
	elif aligner == "bwa":
		if args.r2:
			cmd = f"bwa mem -t {args.thread} {genome_fa} {args.r1} {args.r2} > {sam_out}"
		else:
			cmd = f"bwa mem -t {args.thread} {genome_fa} {args.r1} > {sam_out}"
	elif aligner == "minimap2":
		mmi = os.path.join(aligner_dir, "ref.mmi")
		if args.r2:
			cmd = f"minimap2 -t {args.thread} -ax sr {mmi} {args.r1} {args.r2} > {sam_out}"
		else:
			cmd = f"minimap2 -t {args.thread} -ax sr {mmi} {args.r1} > {sam_out}"
	else:
		print(f"[skip] aligner not supported: {aligner}")
		continue

	run_cmd(cmd)

	print(f"\n[xfq] Filtering reads for {aligner} referencing SAM")

	if aligner == "star":
		sam_out = os.path.join(aligner_dir, "subset-Aligned.out.sam")

	cmd = [
		"python3", "subset.py", "xfq",
		"-g", genome_fa,
		"-s", sam_out,
		"-p", "1",
		"--r1", args.r1
	]

	if args.r2:
		cmd += ["--r2", args.r2]

	cmd += ["-v"]

	run_cmd(" ".join(cmd))

	r1_base = os.path.basename(args.r1)
	r1_root = os.path.splitext(r1_base)[0]
	r1_shrunk = os.path.join(os.path.dirname(args.r1), f"{r1_root}.shrunk.fastq")
	r1_renamed = os.path.join(aligner_dir, f"{r1_root}.shrunk.{aligner}.fastq")
	toolbox.mv(r1_shrunk, r1_renamed)
	if args.r2:
		r2_base = os.path.basename(args.r2)
		r2_root = os.path.splitext(r2_base)[0]
		r2_shrunk = os.path.join(os.path.dirname(args.r2), f"{r2_root}.shrunk.fastq")
		r2_renamed = os.path.join(aligner_dir, f"{r2_root}.shrunk.{aligner}.fastq")
		toolbox.mv(r2_shrunk, r2_renamed)


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
