#!/usr/bin/env python3

import argparse
import gzip
import os
import re
import sys


#########
# tools #
#########

def human_readable_size(num):
	"""Convert a file size in bytes to a human-readable format"""
	suffix = "B"
	for unit in ["", "K", "M", "G", "T", "P", "E", "Z"]:
		if abs(num) < 1024.0:
			return f"{num:3.1f}{unit}{suffix}"
		num /= 1024.0
	return f"{num:.1f}Y{suffix}"


def smart_open_read(filename):
	"""Open a file for reading, support gzipped files"""
	if filename.endswith('.gz'):
		return gzip.open(filename, "rt")
	else:
		return open(filename, "r")
	

def smart_open_write(filename, use_gzip):
	"""Open a file for writing, using gzip if required"""
	if use_gzip:
		return gzip.open(filename, "wt")
	else:
		return open(filename, "w")


#############
# functions #
#############

def write_to_new_genome(cur_header, cur_seqs, scale, outfile, verbose=False):
	"""Helper function to process and write one chromosome's sequence"""
	seq = "".join(cur_seqs)
	new_length = int(len(seq) * scale)
	outfile.write(f">{cur_header}\n")
	outfile.write(seq[:new_length] + "\n")

	if verbose:
		print(f"Processed chromosome {cur_header}:")
		print(f"\tOld length {len(seq)}")
		print(f"\tNew length {new_length}")
	return new_length


def shrink_genome(fin, scale, fout, verbose=False):
	"""
	Reads a FASTA genome file, writes a new FASTA file with first X% of each sequence
	Returns a dictionary mapping chromosome names to new seq lengths
	"""
	new_genome_lengths = {}
	with open(fin, "r") as infile, open(fout, "w") as outfile:
		cur_header = None
		cur_seqs = []

		for line in infile:
			if not line.startswith(">"):
				cur_seqs.append(line.strip())
				continue

			if cur_header is not None:
				new_genome_lengths[cur_header] = write_to_new_genome(cur_header, cur_seqs, scale, outfile, verbose)

			cur_header = line[1:].strip().split()[0]
			cur_seqs = []

		if cur_header is not None:
			new_genome_lengths[cur_header] = write_to_new_genome(cur_header, cur_seqs, scale, outfile, verbose)
	return new_genome_lengths


def get_ref_consumed(cigar):
	"""
	Takes a CIGAR string, returns the reference length elapsed
	Sums the values for operations that consume the reference:
	M, D, N, =, X
	"""
	if cigar == "*" or cigar == "":
		return 0
	ref_len = 0
	# Use regex to find all numbers and their associated operation.
	for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
		if op in "MDN=X":
			ref_len += int(length)
	return ref_len


def validate_reads_from_sam(is_paired, sam):
	"""
	Takes a SAM file
	Returns valid read names still in range of new genome length
	"""
	if is_paired:
		valid_reads = {}
	else:
		valid_reads = set()

	with open(sam, "r") as sam_file:
		for line in sam_file:
			if line.startswith("@"):
				continue

			cols = line.strip().split()
			if len(cols) < 6:
				continue

			qname = cols[0]
			bflag = int(cols[1])
			rname = cols[2]
			pos   = int(cols[3])
			cigar = cols[5]

			if rname not in new_genome_lengths:
				continue

			ref_len = get_ref_consumed(cigar)
			alignment_end = pos - 1 + ref_len

			if alignment_end <= new_genome_lengths[rname]:
				if is_paired:
					if qname not in valid_reads:
						valid_reads[qname] = {"first": False, "second": False}
					if bflag & 64:
						valid_reads[qname]["first"] = True
					elif bflag & 128:
						valid_reads[qname]["second"] = True
					else:
						valid_reads[qname]["first"] = True
						valid_reads[qname]["second"] = True
				else:
					valid_reads.add(qname)
	return valid_reads


def fastq_reader(fp):
	"""Yields one FASTQ record 4 lines at a time"""
	while True:
		try:
			yield [next(fp) for _ in range(4)]
		except StopIteration:
			break


############
# argparse #
############

parser = argparse.ArgumentParser(
	description="Shrink a genome and subset FASTQ files to a given scale")
parser.add_argument("-g", "--genome", type=str, required=True,
	help="Genome FASTA file")
parser.add_argument("-s", "--sam", type=str,
	help="SAM file containing reads alignments")
parser.add_argument("-p", "--scale", type=float, required=True,
	help="Scale to shrink the genome to compared to original")
parser.add_argument("--r1", type=str,
	help="FASTQ file for reads")
parser.add_argument("--r2", type=str, 
	help="FASTQ file for paired-end read2")
parser.add_argument("-z", "--gzip", action="store_true",
	help="Write output FASTQ file(s) in gzipped format")
parser.add_argument("-v", "--verbose", action="store_true",
	help="More verbose output")
args = parser.parse_args()


#########
# check #
#########

if args.scale >= 1:
	sys.exit(f"Error: The scale value must be < 1 (provided scale: {args.scale})")

if args.r2 is not None and args.r1 is None:
	sys.exit("Error: --r2 is provided without --r1")

if args.sam is not None and args.r1 is None:
	sys.exit("Error: SAM file is provided without FASTQ")

if args.r1 is not None and args.sam is None:
	sys.exit("Error: FASTQ file(s) is provided without a SAM file")


##############
# new genome #
##############

if args.verbose:
	print("Starting genome and reads subsetter")
	print(f"\nProcessing genome file from {args.genome}")

shrunken_genome_file = os.path.splitext(args.genome)[0] + ".shrunk.fa"

new_genome_lengths = shrink_genome(args.genome, args.scale, shrunken_genome_file, args.verbose)

if args.verbose:
	print(f"\nShrunken genome written to {shrunken_genome_file}")
	try:
		old_size = os.path.getsize(args.genome)
		new_size = os.path.getsize(shrunken_genome_file)
		print(f"\tOld genome size: {human_readable_size(old_size)}")
		print(f"\tNew genome size: {human_readable_size(new_size)}")
	except Exception as e:
		print("Warning: Unable to retrieve genome file sizes")

if args.sam is None and args.r1 is None and args.r2 is None:
	sys.exit(0)


############
# read SAM #
############

is_paired = args.r2 is not None

if args.verbose:
	print("\nProcessing SAM file for reads subsetting")

valid_reads = validate_reads_from_sam(is_paired, args.sam)

if args.verbose:
	total_valid = len(valid_reads)
	print(f"Unique read names with valid alignments: {total_valid}")


#############
# new FASTQ #
#############

if is_paired:
	out_r1 = os.path.splitext(args.r1)[0] + ".shrunk.fastq"
	out_r2 = os.path.splitext(args.r2)[0] + ".shrunk.fastq"
	if args.gzip:
		out_r1 += ".gz"
		out_r2 += ".gz"
else:
	out_r1 = os.path.splitext(args.r1)[0] + ".shrunk.fastq"
	if args.gzip:
		out_r1 += ".gz"

total_reads = 0
kept_reads = 0

if args.verbose:
	print("\nFiltering reads in FASTQ file(s)")

if is_paired:
	in1  = smart_open_read(args.r1)
	in2  = smart_open_read(args.r2)
	out1 = smart_open_write(out_r1, args.gzip)
	out2 = smart_open_write(out_r2, args.gzip)

	for r1_read, r2_read in zip(fastq_reader(in1), fastq_reader(in2)):
		total_reads += 1
		q1 = r1_read[0].strip()[1:].split()[0]
		q2 = r2_read[0].strip()[1:].split()[0]

		if q1 != q2:
			if args.verbose:
				print(f"Warning: Read names do not match: {q1} vs {q2}. Skipping")
			continue

		if q1 in valid_reads and valid_reads[q1].get("first") and valid_reads[q1].get("second"):
			kept_reads += 1
			out1.write("".join(r1_read))
			out2.write("".join(r2_read))

	in1.close()
	in2.close()
	out1.close()
	out2.close()
else:
	in1  = smart_open_read(args.r1)
	out1 = smart_open_write(out_r1, args.gzip)

	for read in fastq_reader(in1):
		total_reads += 1
		q = read[0].strip()[1:].split()[0]
		if q in valid_reads:
			kept_reads += 1
			out1.write("".join(read))

	in1.close()
	out1.close()

if args.verbose:
	if is_paired:
		print(f"\nTotal read pairs processed from FASTQ: {total_reads}")
		print(f"Read pairs kept: {kept_reads}")
	else:
		print(f"\nTotal reads processed from FASTQ: {total_reads}")
		print(f"Reads kept: {kept_reads}")

try:
	old_r1_size = os.path.getsize(args.r1)
	new_r1_size = os.path.getsize(out_r1)
	print(f"\nFASTQ file: {args.r1}")
	print(f"\tOld size: {human_readable_size(old_r1_size)}")
	print(f"\tNew size: {human_readable_size(new_r1_size)}")
	if is_paired:
		old_r2_size = os.path.getsize(args.r2)
		new_r2_size = os.path.getsize(out_r2)
		print(f"\nFASTQ file: {args.r2}")
		print(f"\tOld size: {human_readable_size(old_r2_size)}")
		print(f"\tNew size: {human_readable_size(new_r2_size)}")
except Exception as e:
	if args.verbose:
		print("Warning: Unable to determine FASTQ file sizes")

print(f"\nGenome and FASTQ subset to {args.scale * 100}%")
