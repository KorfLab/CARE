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
	with open(fin, "r") as infile, open(fout, "w") if fout is not None else open(os.devnull, "w") as outfile:
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
	Takes a CIGAR string, returns the reference length consumed
	Sums the values for operations that consume the reference:
	M, D, N, =, X
	"""
	if cigar == "*" or cigar == "":
		return 0
	ref_len = 0

	for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
		if op in "MDN=X":
			ref_len += int(length)
	return ref_len


def validate_reads_from_sam(is_paired, sam, new_genome_lengths):
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

parser = argparse.ArgumentParser(description="Shrink genome or filter reads with alignment info")
subparsers = parser.add_subparsers(dest="command")

xfa = subparsers.add_parser("xfa",
	help="Shrink a genome FASTA file by scale")
xfa.add_argument("-g", "--genome", required=True,
	help="Genome FASTA file")
xfa.add_argument("-p", "--scale", type=float, required=True,
	help="Scale to shrink compared to original (0 < scale <= 1)")
xfa.add_argument("-v", "--verbose", action="store_true",
	help="More verbose output")

xfq = subparsers.add_parser("xfq",
	help="Filter FASTQ using SAM and genome")
xfq.add_argument("-g", "--genome", required=True,
	help="Genome FASTA file")
xfq.add_argument("-s", "--sam", required=True,
	help="SAM file containing alignment information")
xfq.add_argument("-p", "--scale", type=float, required=True,
	help="Scale to shrink compared to original (0 < scale <= 1)")
xfq.add_argument("--r1", required=True,
	help="FASTQ file for reads")
xfq.add_argument("--r2",
	help="FASTQ file for paired-end read 2")
xfq.add_argument("-z", "--gzip", action="store_true",
	help="Gzip output FASTQ")
xfq.add_argument("-v", "--verbose", action="store_true",
	help="More verbose output")

args = parser.parse_args()

if args.command is None:
	parser.print_help()
	print("\nerror: no subcommand specified, choose one of: {xfa, xfq}")
	sys.exit(1)

if args.scale > 1:
	sys.exit("error: scale must be <= 1")


########
# main #
########

if args.command == "xfa":
	if int(args.scale) == 1:
		sys.exit("Retained original genome")

	if args.verbose:
		print("Starting genome and reads subsetter")
		print(f"\nProcessing genome file from {args.genome}")

	out_genome = os.path.splitext(args.genome)[0] + ".shrunk.fa"

	shrink_genome(args.genome, args.scale, out_genome, args.verbose)

	if args.verbose:
		print(f"\nGenome written to {out_genome}")
	try:
		old_size = os.path.getsize(args.genome)
		new_size = os.path.getsize(out_genome)
		print(f"\tOld genome size: {human_readable_size(old_size)}")
		print(f"\tNew genome size: {human_readable_size(new_size)}")
	except Exception as e:
		print("warning: unable to retrieve genome file sizes")

	print(f"\nGenome subset to {args.scale * 100}%")

elif args.command == "xfq":
	is_paired = args.r2 is not None
	in_genome = args.genome

	if int(args.scale) == 1:
		genome_lengths = shrink_genome(args.genome, 1.0, None, False)
		if args.verbose:
			print("Scale is set to 1, retaining all genome")
	else:
		if args.verbose:
			print("Starting genome and reads subsetter")
			print(f"\nProcessing genome file from {args.genome}")

		out_genome = os.path.splitext(in_genome)[0] + ".shrunk.fa"
		genome_lengths = shrink_genome(in_genome, args.scale, out_genome, args.verbose)

		if args.verbose:
			print(f"\nGenome written to {out_genome}")
		try:
			old_size = os.path.getsize(args.genome)
			new_size = os.path.getsize(out_genome)
			print(f"\tOld genome size: {human_readable_size(old_size)}")
			print(f"\tNew genome size: {human_readable_size(new_size)}")
		except Exception as e:
			print("warning: unable to retrieve genome file sizes")

		print(f"\nGenome subset to {args.scale * 100}%")

	if args.verbose:
		print("\nProcessing SAM file for reads subsetting")

	valid_reads = validate_reads_from_sam(is_paired, args.sam, genome_lengths)

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
		in1 = smart_open_read(args.r1)
		in2 = smart_open_read(args.r2)
		out1 = smart_open_write(out_r1, args.gzip)
		out2 = smart_open_write(out_r2, args.gzip)

		for r1_read, r2_read in zip(fastq_reader(in1), fastq_reader(in2)):
			total_reads += 1
			q1 = r1_read[0].strip()[1:].split()[0]
			q2 = r2_read[0].strip()[1:].split()[0]

			if q1 != q2:
				if args.verbose:
					print(f"warning: read names do not match: {q1} vs {q2}. skipping")
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
		in1 = smart_open_read(args.r1)
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
			print("warning: unable to determine FASTQ file sizes")

	print(f"\nGenome and FASTQ subset to {args.scale * 100}%")
