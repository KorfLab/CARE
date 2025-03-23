#!/usr/bin/env python3

import random
import argparse
import gzip

###############
# Nucleotides #
###############

NUCS = ['A', 'T', 'C', 'G']

###########
# Genomes #
###########

def generate_genome(genome_size, seed):
	random.seed(seed)
	return ''.join(random.choices(NUCS, k=genome_size))


def reverse_complement(seq):
	complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	return ''.join(complement[nuc] for nuc in reversed(seq))

#########
# Reads #
#########

def generate_reads(genome, num_reads, read_length, paired, both_strands, seed):
	random.seed(seed)
	genome_size = len(genome)
	reads = []
	for _ in range(num_reads):
		start = random.randint(0, genome_size - read_length)
		seq = genome[start:start+read_length]
		strand = '+'
		if both_strands and random.choice([True, False]):
			seq = reverse_complement(seq)
			strand = '-'
		qual = 'I' * len(seq)

		header = f"synthetic_read{strand}"

		if paired:
			reads.append((header + '/1', seq, qual))
			start2 = random.randint(0, genome_size - read_length)
			seq2 = genome[start2:start2+read_length]
			if both_strands and random.choice([True, False]):
				seq2 = reverse_complement(seq2)
			reads.append((header + '/2', seq2, qual))
		else:
			reads.append((header, seq, qual))
	return reads

#########
# FastQ #
#########

def write_fastq(filename, reads, gzip_output):
	open_func = gzip.open if gzip_output else open
	if gzip_output:
		filename = filename+".gz"
	mode = 'wt'
	with open_func(filename, mode) as f:
		for i, (header, seq, qual) in enumerate(reads):
			f.write(f"@{header}_{i}\n{seq}\n+\n{qual}\n")

############
# argparse #
############

parser = argparse.ArgumentParser(
	description='Synthetic Genome & Reads Generator')
parser.add_argument('-g', '--genome_size', type=int, default=1000000,
                    help='Genome size (bp)')
parser.add_argument('-n', '--num_reads', type=int, default=100000,
                    help='Number of reads to generate')
parser.add_argument('-l', '--read_length', type=int, default=100,
                    help='Length of each read')
parser.add_argument('-p', '--paired', action='store_true',
                    help='Generate paired-end reads')
parser.add_argument('-b', '--both_strands', action='store_true',
                    help='Generate reads from both strands')
parser.add_argument('-z', '--gzip_output', action='store_true',
                    help='Compress output FASTQ files with gzip')
parser.add_argument('-s', '--seed', type=int, default=1,
                    help='Random seed for reproducibility')
parser.add_argument('-o', '--genome_output', type=str, default='synthetic_genome.fa',
                    help='Output filename for genome')
parser.add_argument('-r', '--reads_output', type=str, default='synthetic_reads.fastq',
                    help='Output filename for reads')
parser.add_argument('-v', '--verbose', action='store_true', default=False,
                    help='Show progress texts')

args = parser.parse_args()

########
# main #
########

if args.verbose:
	print(
		f"Generating {args.genome_size}bp genome and {args.num_reads} of {args.read_length}bp reads")

genome = generate_genome(args.genome_size, args.seed)
with open(args.genome_output, 'w') as gf:
	gf.write(f">synthetic_genome\n")
	for i in range(0, len(genome), 60):
			gf.write(genome[i:i+60] + "\n")

if args.verbose:
	print(f"Saved genome to {args.genome_output}")

reads = generate_reads(genome, args.num_reads, args.read_length,
					   args.paired, args.both_strands, args.seed)
write_fastq(args.reads_output, reads, args.gzip_output)

if args.verbose:
	filename = args.reads_output + ".gz" if args.gzip_output else args.reads_output
	print(f"Saved reads to {filename}")
