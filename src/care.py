#!/usr/bin/env python3

import argparse
import os

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
exp1.add_argument("-k", "--keep", nargs="+", default=["log"],
	help="File extensions to keep after cleanup")

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
exp2.add_argument("-k", "--keep", nargs="+", default=["log"],
	help="File extensions to keep after cleanup")

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
exp3.add_argument("-k", "--keep", nargs="+", default=["log"],
	help="File extensions to keep after cleanup")

# Exp 4: varying read counts
exp4 = subparsers.add_parser('var-rc', help="Align using reads with various read counts")
exp4.add_argument('--r1', nargs='+', required=True,
	help="FASTQ files for varying numbers of reads")
exp4.add_argument('--r2', nargs='+',
	help="FASTQ files for paired-end read 2")
exp4.add_argument('-g', '--genome', required=True,
	help="FASTA file for genome")
exp4.add_argument('-t', '--thread', type=int, default=4,
	help="Number of threads")
exp4.add_argument('-a', '--aligner', nargs='+', required=True,
	help="Aligners to run")
exp4.add_argument('-o', '--output', required=True,
	help="Output directory of experiment")
exp4.add_argument("-k", "--keep", nargs="+", default=["log"],
	help="File extensions to keep after cleanup")

args = parser.parse_args()


########
# main #
########

if args.experiment == "var-gn":
	for aligner in args.aligner:
		aligner_dir = os.path.join(args.output, aligner)
		for genome in args.genome:
			toolbox.sc_fastq(args.r1, args.r2)
			aln = ALIGNER_CLASS_MAP[aligner](genome, aligner_dir, args.r1, args.r2, args.thread, time=True, indexed=False)
			aln.index()
			aln.align()
			aln.cleanup(args.keep)
elif args.experiment == "var-rl":
	for aligner in args.aligner:
		aligner_dir = os.path.join(args.output, aligner)
		indexed = False
		for i, r1 in enumerate(args.r1):
			r2 = args.r2[i] if args.r2 else None
			toolbox.sc_fastq(r1, r2)
			aln = ALIGNER_CLASS_MAP[aligner](args.genome, aligner_dir, r1, r2, args.thread, time=True, indexed=indexed)
			if not indexed:
				aln.index()
				indexed = True
			aln.align()
		aln.cleanup(args.keep)
elif args.experiment == "var-cpu":
	for aligner in args.aligner:
		aligner_dir = os.path.join(args.output, aligner)
		indexed = False
		for thread in args.thread:
			toolbox.sc_fastq(args.r1, args.r2)
			aln = ALIGNER_CLASS_MAP[aligner](args.genome, aligner_dir, args.r1, args.r2, thread, time=True, indexed=indexed)
			if not indexed:
				aln.index()
				indexed = True
			aln.align()
		aln.cleanup(args.keep)
elif args.experiment == "var-rc":
	for aligner in args.aligner:
		aligner_dir = os.path.join(args.output, aligner)
		indexed = False
		for i, r1 in enumerate(args.r1):
			r2 = args.r2[i] if args.r2 else None
			toolbox.sc_fastq(r1, r2)
			aln = ALIGNER_CLASS_MAP[aligner](args.genome, aligner_dir, r1, r2, args.thread, time=True, indexed=indexed)
			if not indexed:
				aln.index()
				indexed = True
			aln.align()
		aln.cleanup(args.keep)
else:
	parser.print_help()
