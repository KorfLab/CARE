# Aligner class

import os
from toolbox import run, cp, only_keep_ext


class Aligner:
	def __init__(self, name, genome, outdir, r1, r2=None, threads=4, time=False):
		self.name = name
		self.genome = genome
		self.outdir = outdir
		self.r1 = r1
		self.r2 = r2
		self.threads = threads

		self.time = time
		self.timelog = os.path.join(outdir, f"{name}.time.log") if time else None

		self.index_dir = os.path.join(outdir, f"{self.name}_index")
		self.ref = os.path.join(self.index_dir, "ref.fa")
		os.makedirs(self.index_dir, exist_ok=True)
		cp(self.genome, self.ref)

		self.sam = None

	def index(self):
		raise NotImplementedError

	def align(self):
		raise NotImplementedError

	def get_sam_path(self):
		return self.sam
	
	def get_timelog_path(self):
		return self.timelog

	def cleanup(self, kept_ext):
		only_keep_ext(self.outdir, kept_ext)


class StarAligner(Aligner):
	def __init__(self, genome, outdir, r1, r2=None, threads=4, time=False):
		super().__init__("star", genome, outdir, r1, r2, threads, time)

	def index(self):
		tmp_dir = os.path.join(self.index_dir, "tmp")
		os.makedirs(tmp_dir, exist_ok=True)

		cmd = [
			"STAR",
			"--runMode", "genomeGenerate",
			"--runThreadN", str(self.threads),
			"--genomeDir", tmp_dir,
			"--genomeFastaFiles", self.ref
		]
		msg = run(cmd)

		recommended = None
		for line in msg.splitlines():
			if "genomeSAindexNbases" in line and "recommended" in line:
				recommended = line.split()[-1]
				break

		if recommended:
			print(f"[STAR] Using recommended genomeSAindexNbases = {recommended}")
			run(["rm", "-rf", tmp_dir])

			cmd = [
				"STAR",
				"--runMode", "genomeGenerate",
				"--runThreadN", str(self.threads),
				"--genomeDir", self.index_dir,
				"--genomeFastaFiles", self.ref,
				"--genomeSAindexNbases", recommended
			]
			run(cmd)
		else:
			print("[STAR] No recommendation found, using default genomeSAindexNbases")
			run(["bash", "-c", f"mv {tmp_dir}/* {self.index_dir}"])
			run(["rm", "-rf", tmp_dir])
		print("[STAR] Indexing complete")

	def align(self):
		prefix = os.path.join(self.outdir, "star_")

		cmd = [
			"STAR",
			"--runMode", "alignReads",
			"--genomeDir", self.index_dir,
			"--runThreadN", str(self.threads),
			"--outFileNamePrefix", prefix,
			"--outSAMtype", "SAM"
		]

		if self.r2:
			cmd += ["--readFilesIn", self.r1, self.r2]
		else:
			cmd += ["--readFilesIn", self.r1]

		if self.time:
			cmd = ["/usr/bin/time", "-v", "-o", self.timelog] + cmd

		run(cmd)
		self.sam = prefix + "Aligned.out.sam"
		print("[STAR] Alignment complete")


class Hisat2Aligner(Aligner):
	def __init__(self, genome, outdir, r1, r2=None, threads=4, time=False):
		super().__init__("hisat2", genome, outdir, r1, r2, threads, time)

	def index(self):
		cmd = [
			"hisat2-build",
			"-p", str(self.threads),
			self.ref,
			self.ref
		]
		run(cmd)
		print("[hisat2] Indexing complete")

	def align(self):
		self.sam = os.path.join(self.outdir, "hisat2.sam")

		if self.r2:
			cmd = [
				"hisat2",
				"-p", str(self.threads),
				"-x", self.ref,
				"-1", self.r1,
				"-2", self.r2,
				"-S", self.sam
			]
		else:
			cmd = [
				"hisat2",
				"-p", str(self.threads),
				"-x", self.ref,
				"-U", self.r1,
				"-S", self.sam
			]

		if self.time:
			cmd = ["/usr/bin/time", "-v", "-o", self.timelog] + cmd

		run(cmd)
		print("[hisat2] Alignment complete")


class Bowtie2Aligner(Aligner):
	def __init__(self, genome, outdir, r1, r2=None, threads=4, time=False):
		super().__init__("bowtie2", genome, outdir, r1, r2, threads, time)

	def index(self):
		cmd = [
			"bowtie2-build",
			"--threads", str(self.threads),
			self.ref,
			self.ref
		]
		run(cmd)
		print("[bowtie2] Indexing complete")

	def align(self):
		self.sam = os.path.join(self.outdir, "bowtie2.sam")

		if self.r2:
			cmd = [
				"bowtie2",
				"-p", str(self.threads),
				"-x", self.ref,
				"-1", self.r1,
				"-2", self.r2,
				"-S", self.sam
			]
		else:
			cmd = [
				"bowtie2",
				"-p", str(self.threads),
				"-x", self.ref,
				"-U", self.r1,
				"-S", self.sam
			]

		if self.time:
			cmd = ["/usr/bin/time", "-v", "-o", self.timelog] + cmd

		run(cmd)
		print("[bowtie2] Alignment complete")


class BwaAligner(Aligner):
	def __init__(self, genome, outdir, r1, r2=None, threads=4, time=False):
		super().__init__("bwa", genome, outdir, r1, r2, threads, time)

	def index(self):
		cmd = [
			"bwa",
			"index",
			self.ref
		]
		run(cmd)
		print("[bwa] Indexing complete")

	def align(self):
		self.sam = os.path.join(self.outdir, "bwa.sam")

		cmd = [
			"bwa",
			"mem",
			"-t",
			str(self.threads),
			self.ref,
			self.r1
		]

		if self.r2:
			cmd.append(self.r2)
		cmd += [">", self.sam]

		if self.time:
			cmd = ["/usr/bin/time", "-v", "-o", self.timelog] + cmd

		run(cmd)
		print("[bwa] Alignment complete")


class Minimap2Aligner(Aligner):
	def __init__(self, genome, outdir, r1, r2=None, threads=4, time=False):
		super().__init__("minimap2", genome, outdir, r1, r2, threads, time)
		self.mmi = os.path.join(self.index_dir, "ref.mmi")

	def index(self):
		cmd = [
			"minimap2",
			"-d", self.mmi,
			self.ref
		]
		run(cmd)
		print("[minimap2] Indexing complete")

	def align(self):
		self.sam = os.path.join(self.outdir, "minimap2.sam")

		cmd = [
			"minimap2",
			"-t", str(self.threads),
			"-ax", "sr",
			self.mmi,
			self.r1
		]

		if self.r2:
			cmd.append(self.r2)
		cmd += [">", self.sam]

		if self.time:
			cmd = ["/usr/bin/time", "-v", "-o", self.timelog] + cmd

		run(cmd)
		print("[minimap2] Alignment complete")
