#!/usr/bin/env python3

import argparse
import os
import sys


#############
# functions #
#############

def run_cmd(cmd):
	"""Run a command and exit when failed"""
	print(f"[run cmd] {cmd}")
	if os.system(cmd) != 0:
		sys.exit(f'FAILED: {cmd}')


def run_cmd_capture(cmd):
	"""Run a command and catch its stdout"""
	print(f"[capture] {cmd}")
	return os.popen(cmd).read()


def index_genome(aligner, genome, thread, aligner_dir):
	"""Index genome inside seperate aligner dir"""
	os.makedirs(aligner_dir, exist_ok=True)
	genome_copy = os.path.join(aligner_dir, "ref.fa")
	run_cmd(f"cp {genome} {genome_copy}")

	if aligner == "star":
		tmp_dir = os.path.join(aligner_dir, "tmp")
		os.makedirs(tmp_dir, exist_ok=True)

		msg = run_cmd_capture(
			f"STAR --runMode genomeGenerate --runThreadN {thread} --genomeDir {tmp_dir} --genomeFastaFiles {genome_copy} 2>&1")
		recommended = None
		for line in msg.splitlines():
			if "genomeSAindexNbases" in line and "recommended" in line:
				recommended = line.split()[-1]
				break

		if recommended is not None:
			print(f"[STAR] Using recommended genomeSAindexNbases = {recommended}")
			run_cmd(f"rm -rf {tmp_dir}")
			run_cmd(
				f"STAR --runMode genomeGenerate --runThreadN {thread} --genomeDir {aligner_dir} --genomeFastaFiles {genome_copy} --genomeSAindexNbases {recommended}")
		else:
			print("[STAR] No recommendation found, using default genomeSAindexNbases")
			run_cmd(f"mv {tmp_dir}/* {aligner_dir}")
			run_cmd(f"rm -rf {tmp_dir}")
	elif aligner == "hisat2":
		run_cmd(f"hisat2-build -p {thread} {genome_copy} {aligner_dir}")
	elif aligner == "bowtie2":
		run_cmd(f"bowtie2-build --threads {thread} {genome_copy} {aligner_dir}")
	elif aligner == "bwa":
		run_cmd(f"bwa index -t {thread} {genome_copy}")
	elif aligner == "minimap2":
		run_cmd(f"minimap2 -d {os.path.join(aligner_dir, 'ref.mmi')} {genome_copy}")
	else:
		print(f"[index_genome] Indexing not supported or needed for {aligner}")


def time_align(aligner, r1, r2, thread, aligner_dir, label):
	"""Perform alignment using aligner in aligner dir"""
	sam = os.path.join(aligner_dir, f"{label}.sam")
	log = os.path.join(aligner_dir, f"{label}.time.txt")
	genome_copy = os.path.join(aligner_dir, "ref.fa")

	if aligner == "star":
		if r2 is None:
			cmd = f"STAR --runMode alignReads --genomeDir {aligner_dir} --readFilesIn {r1} --runThreadN {thread} --outFileNamePrefix {label} --outSAMtype SAM"
		else:
			cmd = f"STAR --runMode alignReads --genomeDir {aligner_dir} --readFilesIn {r1} {r2} --runThreadN {thread} --outFileNamePrefix {label} --outSAMtype SAM"
	elif aligner == "hisat2":
		if r2 is None:
			cmd = f"hisat2 -p {thread} -x {genome_copy} -U {r1} -S {sam}"
		else:
			cmd = f"hisat2 -p {thread} -x {genome_copy} -1 {r1} -2 {r2} -S {sam}"
	elif aligner == "bowtie2":
		if r2 is None:
			cmd = f"bowtie2 -p {thread} -x {genome_copy} -U {r1} -S {sam}"
		else:
			cmd = f"bowtie2 -p {thread} -x {genome_copy} -1 {r1} -2 {r2} -S {sam}"
	elif aligner == "bwa":
		if r2 is None:
			cmd = f"bwa mem -t {thread} {genome_copy} {r1} > {sam}"
		else:
			cmd = f"bwa mem -t {thread} {genome_copy} {r1} {r2} > {sam}"
	elif aligner == "minimap2":
		if r2 is None:
			cmd = f"minimap2 -t {thread} -ax sr {os.path.join(aligner_dir, 'ref.mmi')} {r1} > {sam}"
		else:
			cmd = f"minimap2 -t {thread} -ax sr {os.path.join(aligner_dir, 'ref.mmi')} {r1} {r2} > {sam}"
	else:
		print(f"Unsupported aligner: {aligner}")
		return

	full_cmd = f"/usr/bin/time -v {cmd} 2> {log}"
	run_cmd(full_cmd)


############
# argparse #
############

parser = argparse.ArgumentParser(
	description="Run alignment experiments with various inputs")
subparsers = parser.add_subparsers(dest='experiment', help='Experiment type')

# Exp 1: varying genome sizs
exp1 = subparsers.add_parser('var-gn', help="Align against genomes with various sizes")
exp1.add_argument('--r1', required=True,
                  help="FASTQ file for reads")
exp1.add_argument('--r2',
                  help="FASTQ file for paired-end read 2")
exp1.add_argument('-g', '--genome', nargs='+', required=True,
                  help="FASTA file for genome of varying sizes")
exp1.add_argument('-t', '--thread', type=int, default=4,
                  help="Number of threads")
exp1.add_argument('-a', '--aligner', nargs='+', required=True,
                  help="Aligners to run")
exp1.add_argument('-o', '--output', required=True,
                  help="Output directory of experiment")

# Exp 2: varying read lengths
exp2 = subparsers.add_parser('var-rl', help="Align using reads with various read lengths")
exp2.add_argument('--r1', nargs='+', required=True,
                  help="FASTQ files for reads of varying read lengths")
exp2.add_argument('--r2', nargs='+',
                  help="FASTQ files for paired-end read 2")
exp2.add_argument('-g', '--genome', required=True,
                  help="FASTA file for genome")
exp2.add_argument('-t', '--thread', type=int, default=4,
                  help="Number of threads")
exp2.add_argument('-a', '--aligner', nargs='+', required=True,
                  help="Aligners to run")
exp2.add_argument('-o', '--output', required=True,
                  help="Output directory of experiment")

# Exp 3: varying number of threads
exp3 = subparsers.add_parser('var-cpu', help="Align using various numbers of threads")
exp3.add_argument('--r1', required=True,
                  help="FASTQ file for reads")
exp3.add_argument('--r2',
                  help="FASTQ file for paired-end read 2")
exp3.add_argument('-g', '--genome', required=True,
                  help="FASTA file for genome")
exp3.add_argument('-t', '--thread', type=int, nargs='+', required=True,
                  help="Various number of threads")
exp3.add_argument('-a', '--aligner', nargs='+', required=True,
                  help="Aligners to run")
exp3.add_argument('-o', '--output', required=True,
                  help="Output directory of experiment")

# Exp 4: varying coverage
exp4 = subparsers.add_parser('var-cov', help="Align using reads with various coverages")
exp4.add_argument('--r1', nargs='+', required=True,
                  help="FASTQ files for reads of varying coverage")
exp4.add_argument('--r2',
                  help="FASTQ files for paired-end read 2")
exp4.add_argument('-g', '--genome', required=True,
                  help="FASTA file for genome")
exp4.add_argument('-t', '--thread', type=int, default=4,
                  help="Number of threads")
exp4.add_argument('-a', '--aligner', nargs='+', required=True,
                  help="Aligners to run")
exp4.add_argument('-o', '--output', required=True,
                  help="Output directory of experiment")

args = parser.parse_args()


########
# main #
########

for aligner in args.aligner:
	aligner_dir = os.path.join(args.output, aligner)
	if args.experiment == "var-gn":
		for genome in args.genome:
			index_genome(aligner, genome, args.thread, aligner_dir)
	elif args.experiment == "var-cpu":
		index_genome(aligner, args.genome, 4, aligner_dir)
	else:
		index_genome(aligner, args.genome, args.thread, aligner_dir)

if args.experiment == "var-gn":
	for aligner in args.aligner:
		aligner_dir = os.path.join(args.output, aligner)
		for genome in args.genome:
			label = f"{os.path.basename(genome)}_{aligner}"
			time_align(aligner, args.r1, args.r2, args.thread, aligner_dir, label)
elif args.experiment == "var-rl":
	for aligner in args.aligner:
		aligner_dir = os.path.join(args.output, aligner)
		if args.r2 is None:
			for reads in args.r1:
				label = f"{os.path.basename(reads)}_{aligner}"
				time_align(aligner, reads, None, args.thread, aligner_dir, label)
		else:
			if len(args.r1) != len(args.r2):
				sys.exit(f'ERROR: Unmatched number of pair-ended reads')
			for i in range(len(args.r1)):
				reads1 = args.r1[i]
				reads2 = args.r2[i]
				label = f"{os.path.basename(reads1)}_{aligner}"
				time_align(aligner, reads1, reads2, args.thread, aligner_dir, label)
elif args.experiment == "var-cpu":
	for aligner in args.aligner:
		aligner_dir = os.path.join(args.output, aligner)
		for thread in args.thread:
			label = f"{aligner}_t{thread}"
			time_align(aligner, args.r1, args.r2, thread, aligner_dir, label)
elif args.experiment == "var-cov":
	for aligner in args.aligner:
		aligner_dir = os.path.join(args.output, aligner)
		if args.r2 is None:
			for reads in args.r1:
				label = f"{os.path.basename(reads)}_{aligner}"
				time_align(aligner, reads, None, args.thread, aligner_dir, label)
		else:
			if len(args.r1) != len(args.r2):
				sys.exit(f'ERROR: Unmatched number of pair-ended reads')
			for reads1 in args.r1:
				for reads2 in args.r2:
					label = f"{os.path.basename(reads1)}_{aligner}"
					time_align(aligner, reads1, reads2, args.thread, aligner_dir, label)
else:
	parser.print_help()
