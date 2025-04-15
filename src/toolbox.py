# toolbox.py

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


def only_keep_ext(directory, ext):
	"""Delete everything in a directory except files with the extension"""
	if not ext.startswith("."):
		ext = "." + ext

	at_least_one = False

	for item in os.listdir(directory):
		path = os.path.join(directory, item)

		if os.path.isdir(path):
			shutil.rmtree(path)
			print(f"[only keep ext] Removed directory: {path}")
		elif not item.endswith(ext):
			os.remove(path)
			print(f"[only keep ext] Removed file: {path}")
		else:
			at_least_one = True

	if not at_least_one:
		print(f"[only keep ext] WARNING: No *{ext} files found in {directory}")


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

		if recommended is not None:
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
