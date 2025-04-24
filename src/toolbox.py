# toolbox.py

import datetime
import gzip
import os
import shutil
import subprocess
import sys


def run(cmd):
	"""Runs list-based command and exits when failed, returns stdout"""
	cmd_str = ' '.join(cmd)
	print(f"[run] {cmd_str}")

	if ">" in cmd_str:
		result = subprocess.run(cmd_str, shell=True, capture_output=True, text=True)
	else:
		result = subprocess.run(cmd, capture_output=True, text=True)

	if result.returncode != 0:
		print(f"[run] FAILED: {cmd_str}")
		print(f"[run] ERROR: {result.stderr}")
		sys.exit(1)

	return result.stdout


def cp(src, dst):
	"""Copy file or directory from src to dst"""
	if not os.path.exists(src):
		raise FileNotFoundError(f"[cp] ERROR: Source '{src}' does not exist")
	if os.path.isdir(src):
		if os.path.isdir(dst):
			dst = os.path.join(dst, os.path.basename(src))
		if os.path.exists(dst):
			raise FileExistsError(f"[cp] ERROR: Target directory '{dst}' already exists")
		shutil.copytree(src, dst)
		print(f"[cp] Copied directory: {src} -> {dst}")
		return dst
	else:
		if os.path.isdir(dst):
			dst = os.path.join(dst, os.path.basename(src))
		if os.path.isfile(dst):
			raise FileExistsError(f"[cp] ERROR: Target file '{dst}' already exists")
		result = shutil.copy(src, dst)
		print(f"[cp] Copied file: {src} -> {result}")
		return result


def mv(src, dst):
	"""Move file or directory from src to dst"""
	if not os.path.exists(src):
		raise FileNotFoundError(f"[mv] ERROR: Source '{src}' does not exist")

	if os.path.isdir(dst):
		dst = os.path.join(dst, os.path.basename(src))

	if os.path.exists(dst):
		raise FileExistsError(f"[mv] ERROR: Target '{dst}' already exists")

	result = shutil.move(src, dst)
	print(f"[mv] Moved: {src} -> {result}")
	return result


def only_keep_ext(directory, ext):
	"""Delete everything in a directory except files with the given extensions"""
	if isinstance(ext, str):
		ext = [ext]

	ext = [f".{e}" if not e.startswith(".") else e for e in ext]

	at_least_one = False

	for item in os.listdir(directory):
		path = os.path.join(directory, item)

		if os.path.isdir(path):
			shutil.rmtree(path)
			print(f"[only keep ext] Removed directory: {path}")
		elif not any(item.endswith(e) for e in ext):
			os.remove(path)
			print(f"[only keep ext] Removed file: {path}")
		else:
			at_least_one = True

	if not at_least_one:
		print(f"[only keep ext] WARNING: No matching {' '.join(ext)} files found in {directory}")


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
	

def smart_open_append(filename):
	"""Open a file for appending, support gzipped files"""
	if filename.endswith('.gz'):
		return gzip.open(filename, "at")
	else:
		return open(filename, "a")


def fastq_reader(fp):
	"""Yields one FASTQ record 4 lines at a time"""
	while True:
		try:
			yield [next(fp) for _ in range(4)]
		except StopIteration:
			break


def sc_fastq(file1, file2=None):
	"""Sanity check single or paired FASTQ"""

	with smart_open_read(file1) as f1:
		lines1 = sum(1 for _ in f1)

	if lines1 % 4 != 0:
		print(f"[sc-fastq] ERROR: {file1} has {lines1} lines, which is not a multiple of 4")
		sys.exit(1)

	reads1 = lines1 // 4
	print(f"[sc-fastq] {file1} has {reads1} valid reads")

	if file2:
		with smart_open_read(file2) as f2:
			lines2 = sum(1 for _ in f2)

		if lines2 % 4 != 0:
			print(f"[sc-fastq] ERROR: {file2} has {lines2} lines, which is not a multiple of 4")
			sys.exit(1)

		reads2 = lines2 // 4

		if reads1 != reads2:
			print(f"[sc-fastq] ERROR: Mismatched read counts:\n"
			      f"\t{file1} - {reads1} reads\n"
			      f"\t{file2} - {reads2} reads")
			sys.exit(1)

		print(f"[sc-fastq] {file2} has {reads2} valid reads")
		print(f"[sc-fastq] SC Complete: {file1} and {file2} validated")
	else:
		print(f"[sc-fastq] SC Complete: {file1} validated")


def get_num_reads(fq):
	"""Return the number of reads in a FASTQ file"""
	with smart_open_read(fq) as f:
		return sum(1 for _ in f) // 4


def timestamp():
	"""Return timestamp string of YYYYMMDD-HHMMSS-mmm"""
	return datetime.datetime.now().strftime("%Y%m%d-%H%M%S-%f")[:-3]


def get_genome_length(fasta_file):
	"""Return total number of bases in a FASTA file"""
	total = 0
	with open(fasta_file, "r") as fa:
		for line in fa:
			if not line.startswith(">"):
				total += len(line.strip())
	return total


def get_read_length(fastq_file):
	"""Return the length of the first read in a FASTQ file"""
	with smart_open_read(fastq_file) as fq:
		next(fq)
		seq = next(fq).strip()
		return len(seq)


def read_fasta(fasta_file):
	"""Reads a FASTA file and returns a dict of {header: sequence}"""
	sequences = {}
	header = None
	seq_lines = []

	with smart_open_read(fasta_file) as fa:
		for line in fa:
			if line.startswith(">"):
				if header:
					sequences[header] = "".join(seq_lines)
				header = line[1:].strip().split()[0]
				seq_lines = []
			else:
				seq_lines.append(line.strip())

	if header:
		sequences[header] = "".join(seq_lines)

	return sequences
