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


def fastq_reader(fp):
	"""Yields one FASTQ record 4 lines at a time"""
	while True:
		try:
			yield [next(fp) for _ in range(4)]
		except StopIteration:
			break


def count_lines(fp):
	"""Count lines in fp"""
	count = 0
	for _ in fp:
		count += 1
	return count


def sc_fastq(file1, file2=None):
	"""Sanity check single or paired FASTQ"""

	with smart_open_read(file1) as f1:
		lines1 = count_lines(f1)

	if lines1 % 4 != 0:
		print(f"[sc-fastq] ERROR: {file1} has {lines1} lines, which is not a multiple of 4")
		sys.exit(1)

	reads1 = lines1 // 4
	print(f"[sc-fastq] {file1} has {reads1} valid reads")

	if file2:
		with smart_open_read(file2) as f2:
			lines2 = count_lines(f2)

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


def timestamp():
	"""Return timestamp string of YYYYMMDD-HHMMSS-mmm"""
	return datetime.datetime.now().strftime("%Y%m%d-%H%M%S-%f")[:-3]


def index(aligner, genome, threads, index_dir):
	"""Index a genome for an aligner"""
	os.makedirs(index_dir, exist_ok=True)
	genome_copy = os.path.join(index_dir, "ref.fa")
	shutil.copy(genome, genome_copy)

	if aligner == "star":
		tmp_dir = os.path.join(index_dir, "tmp")
		os.makedirs(tmp_dir, exist_ok=True)

		cmd = [
			"STAR",
			"--runMode", "genomeGenerate",
			"--runThreadN", str(threads),
			"--genomeDir", tmp_dir,
			"--genomeFastaFiles", genome_copy
		]
		msg = run(cmd)

		recommended = None
		for line in msg.splitlines():
			if "genomeSAindexNbases" in line and "recommended" in line:
				recommended = line.split()[-1]
				break

		if recommended:
			print(f"[{aligner}] Using recommended genomeSAindexNbases = {recommended}")
			run(["rm", "-rf", tmp_dir])

			cmd = [
				"STAR",
				"--runMode", "genomeGenerate",
				"--runThreadN", str(threads),
				"--genomeDir", index_dir,
				"--genomeFastaFiles", genome_copy,
				"--genomeSAindexNbases", recommended
			]
			run(cmd)
		else:
			print(f"[{aligner}] No recommendation found, using default genomeSAindexNbases")
			shutil.move(tmp_dir, index_dir)
		print(f"[{aligner}] Indexing complete")
	elif aligner == "hisat2":
		cmd = [
			"hisat2-build",
			"-p", str(threads),
			genome_copy,
			genome_copy
		]
		run(cmd)
		print(f"[{aligner}] Indexing complete")
	elif aligner == "bowtie2":
		cmd = [
			"bowtie2-build",
			"--threads", str(threads),
			genome_copy,
			genome_copy
		]
		run(cmd)
		print(f"[{aligner}] Indexing complete")
	elif aligner == "bwa":
		cmd = [
			"bwa",
			"index",
			genome_copy
		]
		run(cmd)
		print(f"[{aligner}] Indexing complete")
	elif aligner == "minimap2":
		cmd = [
			"minimap2",
			"-d", os.path.join(index_dir, "ref.mmi"),
			genome_copy
		]
		run(cmd)
		print(f"[{aligner}] Indexing complete")
	# add new aligner command here
	else:
		print(f"[index] Indexing not supported for aligner: {aligner}")


def align(aligner, r1, r2, threads, index_dir, prefix):
	"""Run alignment using an aligner"""
	ref = os.path.join(index_dir, "ref.fa")
	if aligner == "star":
		if r2:
			cmd = [
				"STAR",
				"--runMode", "alignReads",
				"--genomeDir", index_dir,
				"--readFilesIn", r1, r2,
				"--runThreadN", str(threads),
				"--outFileNamePrefix", prefix,
				"--outSAMtype", "SAM"
			]
		else:
			cmd = [
				"STAR",
				"--runMode", "alignReads",
				"--genomeDir", index_dir,
				"--readFilesIn", r1,
				"--runThreadN", str(threads),
				"--outFileNamePrefix", prefix,
				"--outSAMtype", "SAM"
			]
		run(cmd)
		print(f"[{aligner}] Alignment complete")
	elif aligner == "hisat2":
		if r2:
			cmd = [
				"hisat2",
				"-p", str(threads),
				"-x", ref,
				"-1", r1,
				"-2", r2,
				"-S", f"{prefix}.sam"
			]
		else:
			cmd = [
				"hisat2",
				"-p", str(threads),
				"-x", ref,
				"-U", r1,
				"-S", f"{prefix}.sam"
			]
		run(cmd)
		print(f"[{aligner}] Alignment complete")
	elif aligner == "bowtie2":
		if r2:
			cmd = [
				"bowtie2",
				"-p", str(threads),
				"-x", ref,
				"-1", r1,
				"-2", r2,
				"-S", f"{prefix}.sam"
			]
		else:
			cmd = [
				"bowtie2",
				"-p", str(threads),
				"-x", ref,
				"-U", r1,
				"-S", f"{prefix}.sam"
			]
		run(cmd)
		print(f"[{aligner}] Alignment complete")
	elif aligner == "bwa":
		if r2:
			cmd = [
				"bwa",
				"mem",
				"-t", str(threads),
				ref,
				r1,
				r2
			]
			cmd += [f"> {prefix}.sam"]
		else:
			cmd = [
				"bwa",
				"mem",
				"-t", str(threads),
				ref,
				r1
			]
			cmd += [f"> {prefix}.sam"]
		run(cmd)
		print(f"[{aligner}] Alignment complete")
	elif aligner == "minimap2":
		mmi = os.path.join(index_dir, "ref.mmi")
		if r2:
			cmd = [
				"minimap2",
				"-t", str(threads),
				"-ax", "sr",
				mmi,
				r1,
				r2
			]
			cmd += [f"> {prefix}.sam"]
		else:
			cmd = [
				"minimap2",
				"-t", str(threads),
				"-ax", "sr",
				mmi,
				r1
			]
			cmd += [f"> {prefix}.sam"]
		run(cmd)
		print(f"[{aligner}] Alignment complete")
	# add new aligner command here
	else:
		print(f"[align] Alignment not supported for aligner: {aligner}")
