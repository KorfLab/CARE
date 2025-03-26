#!/usr/bin/env python3
import argparse
import os
import gzip
import random
import sys


#########
# Tools #
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


def remove_extensions(filename, extensions):
	"""Remove one of the given extensions from filename"""
	for ext in extensions:
		if filename.endswith(ext):
			return filename[: -len(ext)]
	return filename


def validate_pair_fastq(file1, file2, verbose=False):
	"""
	Ensure that both paired-end FASTQ files have the same number of reads
	Counts the total lines in each file and verifies that each is a multiple of 4,
	and that the number of reads (lines/4) in both files are equal
	"""
	f1 = smart_open_read(file1)
	f2 = smart_open_read(file2)
	count1 = 0
	count2 = 0
	for _ in f1:
		count1 += 1
	for _ in f2:
		count2 += 1
	f1.close()
	f2.close()

	if count1 % 4 != 0 or count2 % 4 != 0:
		sys.exit("Error: One of the input FASTQ files does not contain a multiple of 4 lines")
	reads1 = count1 // 4
	reads2 = count2 // 4
	if reads1 != reads2:
		sys.exit(f"Error: Paired FASTQ files have a different number of reads: {reads1} versus {reads2}")
	if verbose:
		print(f"Paired files validated: {reads1} read pairs found in each file")
	return reads1


#############
# Reservior #
#############

def reservoir_sample_single(fp, k):
	"""
	Reservoir sampling for a single-end FASTQ file
	Each read is 4 consecutive lines
	Returns a list of selected reads and the total read count
	"""
	reservoir = []
	total = 0
	while True:
		try:
			read = [next(fp) for _ in range(4)]
		except StopIteration:
			break
		total += 1
		if total <= k:
			reservoir.append(read)
		else:
			# With probability k/total, replace a read in the reservoir
			r = random.randint(1, total)
			if r <= k:
				reservoir[r - 1] = read
	return reservoir, total


def reservoir_sample_paired(fp1, fp2, k):
	"""
	Reservoir sampling for paired-end FASTQ files
	Each read pair comprises one read from each file
	Returns two reservoirs and the total number of read pairs
	"""
	reservoir1 = []
	reservoir2 = []
	total = 0
	while True:
		try:
			read1 = [next(fp1) for _ in range(4)]
			read2 = [next(fp2) for _ in range(4)]
		except StopIteration:
			break
		total += 1
		if total <= k:
			reservoir1.append(read1)
			reservoir2.append(read2)
		else:
			r = random.randint(1, total)
			if r <= k:
				reservoir1[r - 1] = read1
				reservoir2[r - 1] = read2
	return reservoir1, reservoir2, total


############
# argparse #
############

parser = argparse.ArgumentParser(description="Shrink FASTQ file by randomly subsetting")
parser.add_argument("--r1", required=True,
	help="Path to input FASTQ file for reads 1")
parser.add_argument("--r2", 
	help="Path to input FASTQ file for reads 2 for paired-end")
parser.add_argument("-n", "--numReads", type=int, required=True,
	help="Number of reads to keep in the output")
parser.add_argument("-s", "--seed", type=int, default=1,
	help="Set seed for random number generator")
parser.add_argument("-z", "--gzip", action="store_true",
	help="Output file in gzip-compressed format")
parser.add_argument("--sort", action="store_true",
	help="Sort the output by header. Only for NCBI reads")
parser.add_argument("-v", "--verbose", action="store_true",
	help="More verbose output")
args = parser.parse_args()


######
# in #
######

random.seed(args.seed)

r1_size = os.path.getsize(args.r1)
if args.r2:
	r2_size = os.path.getsize(args.r2)
if args.verbose:
	print("Starting minifq")
	print("Input file(s):")
	print(f"\t{args.r1} - {human_readable_size(r1_size)}")
	if args.r2:
		print(f"\t{args.r2} - {human_readable_size(r2_size)}")

if args.r2:
	total_pairs = validate_pair_fastq(args.r1, args.r2, args.verbose)
	fin1 = smart_open_read(args.r1)
	fin2 = smart_open_read(args.r2)
	reservoir1, reservoir2, total_reads = reservoir_sample_paired(fin1, fin2, args.numReads)
	fin1.close()
	fin2.close()
	if args.verbose:
		print(f"Total read pairs in input: {total_reads}")
		print(f"Selecting {args.numReads} read pairs randomly")
else:
	fin1 = smart_open_read(args.r1)
	reservoir, total_reads = reservoir_sample_single(fin1, args.numReads)
	fin1.close()
	if args.verbose:
		print(f"Total reads in input: {total_reads}")
		print(f"Selecting {args.numReads} reads randomly")


#######
# out #
#######

extensions_to_remove = [".fq", ".fastq", ".fq.gz", ".fastq.gz"]
if args.r2:
	if args.sort:
		if args.verbose:
			print("Sorting new reads by header")
		paired = list(zip(reservoir1, reservoir2))
		paired.sort(key=lambda pair: int(pair[0][0].strip().split()[0].split('.')[1]))
		reservoir1, reservoir2 = zip(*paired)
		reservoir1 = list(reservoir1)
		reservoir2 = list(reservoir2)
	base1 = remove_extensions(os.path.basename(args.r1), extensions_to_remove)
	base2 = remove_extensions(os.path.basename(args.r2), extensions_to_remove)
	out1 = base1 + ".minifq.fastq"
	out2 = base2 + ".minifq.fastq"
	if args.gzip:
		out1 += ".gz"
		out2 += ".gz"
	if args.verbose:
		print(f"Writing output to:\n\t{out1}\n\t{out2}")
	fout1 = smart_open_write(out1, args.gzip)
	fout2 = smart_open_write(out2, args.gzip)
	for read in reservoir1:
		fout1.write("".join(read))
	for read in reservoir2:
		fout2.write("".join(read))
	fout1.close()
	fout2.close()
	out1_size = os.path.getsize(out1)
	out2_size = os.path.getsize(out2)
	if args.verbose:
		print("Output file sizes:")
		print(f"\t{out1} - {human_readable_size(out1_size)}")
		print(f"\t{out2} - {human_readable_size(out2_size)}")
		print("Finished processing paired-end files")
else:
	if args.sort:
		if args.verbose:
			print("Sorting new reads by header")
		reservoir = sorted(reservoir, key=lambda read: int(read[0].strip().split()[0].split('.')[1]))
	base1 = remove_extensions(os.path.basename(args.r1), extensions_to_remove)
	out_file = base1 + ".minifq.fastq"
	if args.gzip:
		out_file += ".gz"
	if args.verbose:
		print(f"Writing output to: {out_file}")
	fout = smart_open_write(out_file, args.gzip)
	for read in reservoir:
		fout.write("".join(read))
	fout.close()
	out_size = os.path.getsize(out_file)
	if args.verbose:
		print(f"Output file {out_file} size: {human_readable_size(out_size)}")
		print("Finished processing single-end file")
