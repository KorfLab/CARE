#!/usr/bin/env python3

import argparse
import os
import random
import toolbox


#########
# funcs #
#########

def remove_extensions(filename, extensions):
	"""Remove one of the given extensions from filename"""
	for ext in extensions:
		if filename.endswith(ext):
			return filename[: -len(ext)]
	return filename


#############
# reservior #
#############

def reservoir_sample_single(fp, k):
	"""
	Reservoir sampling for a single-end FASTQ file
	Each read is 4 consecutive lines
	Returns a list of selected reads and the total read count
	"""
	reservoir = []
	total = 0

	for read in toolbox.fastq_reader(fp):
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

	for read1, read2 in zip(toolbox.fastq_reader(fp1), toolbox.fastq_reader(fp2)):
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
parser.add_argument("-o", "--output", default=".",
	help="Directory to save output FASTQ files (default: current directory)")
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

outdir = args.output
os.makedirs(outdir, exist_ok=True)


######
# in #
######

random.seed(args.seed)

r1_size = os.path.getsize(args.r1)
if args.r2:
	r2_size = os.path.getsize(args.r2)

if args.verbose:
	print("[minifq] Starting minifq")
	print("[minifq] Input file(s):")
	print(f"\t{args.r1} - {toolbox.human_readable_size(r1_size)}")
	if args.r2:
		print(f"\t{args.r2} - {toolbox.human_readable_size(r2_size)}")

if args.r2:
	toolbox.sc_fastq(args.r1, args.r2)
	fin1 = toolbox.smart_open_read(args.r1)
	fin2 = toolbox.smart_open_read(args.r2)
	reservoir1, reservoir2, total_reads = reservoir_sample_paired(fin1, fin2, args.numReads)
	fin1.close()
	fin2.close()

	if args.verbose:
		print(f"[minifq] Total read pairs in input: {total_reads}")
		print(f"[minifq] Selecting {args.numReads} read pairs randomly")
else:
	toolbox.sc_fastq(args.r1)
	fin1 = toolbox.smart_open_read(args.r1)
	reservoir, total_reads = reservoir_sample_single(fin1, args.numReads)
	fin1.close()

	if args.verbose:
		print(f"[minifq] Total reads in input: {total_reads}")
		print(f"[minifq] Selecting {args.numReads} reads randomly")


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
	out1 = os.path.join(outdir, base1 + ".minifq.fastq")
	out2 = os.path.join(outdir, base2 + ".minifq.fastq")

	if args.gzip:
		out1 += ".gz"
		out2 += ".gz"

	if args.verbose:
		print(f"Writing output to:\n\t{out1}\n\t{out2}")

	fout1 = toolbox.smart_open_write(out1, args.gzip)
	fout2 = toolbox.smart_open_write(out2, args.gzip)

	for read in reservoir1:
		fout1.write("".join(read))
	for read in reservoir2:
		fout2.write("".join(read))

	fout1.close()
	fout2.close()

	out1_size = os.path.getsize(out1)
	out2_size = os.path.getsize(out2)

	if args.verbose:
		print("[minifq] Output file sizes:")
		print(f"\t{out1} - {toolbox.human_readable_size(out1_size)}")
		print(f"\t{out2} - {toolbox.human_readable_size(out2_size)}")
		print("[minifq] Finished processing paired-end files")
else:
	if args.sort:
		if args.verbose:
			print("[minifq] Sorting new reads by header")
		reservoir = sorted(reservoir, key=lambda read: int(read[0].strip().split()[0].split('.')[1]))

	base1 = remove_extensions(os.path.basename(args.r1), extensions_to_remove)
	out_file = os.path.join(outdir, base1 + ".minifq.fastq")

	if args.gzip:
		out_file += ".gz"

	if args.verbose:
		print(f"[minifq] Writing output to: {out_file}")

	fout = toolbox.smart_open_write(out_file, args.gzip)

	for read in reservoir:
		fout.write("".join(read))
	fout.close()
	out_size = os.path.getsize(out_file)

	if args.verbose:
		print(f"[minifq] Output file size:")
		print(f"\t{out_file} - {toolbox.human_readable_size(out_size)}")
		print("[minifq] Finished processing single-end file")
