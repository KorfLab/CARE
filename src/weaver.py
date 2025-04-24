#!/usr/bin/env python3

import argparse
import os
import random
import sys

import toolbox

############
# argparse #
############

parser = argparse.ArgumentParser(description="Generate or extend FASTQ files to a set number of reads")
subparsers = parser.add_subparsers(dest="command", help="weaver modes: reuse, extend, synth")

reuse = subparsers.add_parser("reuse",
	help="Reuse existing reads")
reuse.add_argument("--r1", required=True,
	help="FASTQ file for reads")
reuse.add_argument("--r2",
	help="FASTQ file for paired-end read 2")
reuse.add_argument("-n", "--numReads", type=int, required=True,
	help="Target total number of reads")
reuse.add_argument("-o", "--output", default=".",
	help="Directory to save output FASTQ files (default: current directory)")
reuse.add_argument("-s", "--seed", type=int, default=1,
	help="Random seed (default: 1)")
reuse.add_argument("-v", "--verbose", action="store_true",
	help="Verbose logging")

extend = subparsers.add_parser("extend", help="Extend existing reads with synthetic reads from genome")
extend.add_argument("--r1", required=True,
	help="FASTQ file for reads")
extend.add_argument("--r2",
	help="FASTQ file for paired-end read 2")
extend.add_argument("-g", "--genome", required=True,
	help="Genome FASTA file to use for read generation")
extend.add_argument("-n", "--numReads", type=int, required=True,
	help="Target total number of reads")
extend.add_argument("-o", "--output", default=".",
	help="Directory to save output FASTQ files (default: current directory)")
extend.add_argument("-s", "--seed", type=int, default=1,
	help="Random seed (default: 1)")
extend.add_argument("-v", "--verbose", action="store_true",
	help="Verbose logging")

synth = subparsers.add_parser("synth", help="Generate synthetic reads from genome only")
synth.add_argument("-g", "--genome", required=True,
	help="Genome FASTA file to use for read generation")
synth.add_argument("-k", "--readLength", type=int, required=True,
	help="Read length of synthetic reads")
synth.add_argument("-n", "--numReads", type=int, required=True,
	help="Target number of reads")
synth.add_argument("-p", "--paired", action="store_true",
	help="Generate paired-end reads")
synth.add_argument("-o", "--output", default=".",
	help="Directory to save output FASTQ files (default: current directory)")
synth.add_argument("--prefix", default="synthetic",
	help="Prefix for output FASTQ files")
synth.add_argument("-s", "--seed", type=int, default=1,
	help="Random seed (default: 1)")
synth.add_argument("-v", "--verbose", action="store_true",
	help="Verbose logging")

args = parser.parse_args()


########
# main #
########

if args.command == "reuse":
	toolbox.sc_fastq(args.r1, args.r2)

	if args.output:
		os.makedirs(args.output, exist_ok=True)

	random.seed(args.seed)
	target_n = args.numReads

	base_r1 = os.path.splitext(os.path.basename(args.r1))[0]
	r1_out = os.path.join(args.output, f"{base_r1}.weaver.fastq")
	toolbox.cp(args.r1, r1_out)

	if args.r2:
		base_r2 = os.path.splitext(os.path.basename(args.r2))[0]
		r2_out = os.path.join(args.output, f"{base_r2}.weaver.fastq")
		toolbox.cp(args.r2, r2_out)

	if args.verbose:
		print(f"[weaver] Copied original reads to:\n\t{r1_out}" + (f"\n\t{r2_out}" if args.r2 else ""))

	with toolbox.smart_open_read(args.r1) as f:
		orig_reads = sum(1 for _ in f) // 4

	if args.verbose:
		print(f"[weaver] Original reads: {orig_reads}")

	if orig_reads >= target_n:
		print(f"[weaver] Requested reads ({target_n}) already satisfied by original")
		sys.exit(0)

	num_to_add = target_n - orig_reads
	if args.verbose:
		print(f"[weaver] Reusing reads to add {num_to_add} more")

	with toolbox.smart_open_read(args.r1) as fin1:
		r1_pool = [read for read in toolbox.fastq_reader(fin1)]

	if args.r2:
		with toolbox.smart_open_read(args.r2) as fin2:
			r2_pool = [read for read in toolbox.fastq_reader(fin2)]
	
	read_len = toolbox.get_read_length(args.r1)

	with toolbox.smart_open_append(r1_out) as fout1:
		if args.r2:
			with toolbox.smart_open_append(r2_out) as fout2:
				for i in range(num_to_add):
					idx = random.randint(0, len(r1_pool) - 1)
					r1 = r1_pool[idx][:]
					r2 = r2_pool[idx][:]

					hdr = f"WEAVERSYNTH.{i+1} {i+1} length={read_len}"
					r1[0] = f"@{hdr}\n"
					r1[2] = f"+{hdr}\n"
					r2[0] = f"@{hdr}\n"
					r2[2] = f"+{hdr}\n"

					fout1.write("".join(r1))
					fout2.write("".join(r2))
		else:
			for i in range(num_to_add):
				idx = random.randint(0, len(r1_pool) - 1)
				r1 = r1_pool[idx][:]

				hdr = f"WEAVERSYNTH.{i+1} {i+1} length={read_len}"
				r1[0] = f"@{hdr}\n"
				r1[2] = f"+{hdr}\n"

				fout1.write("".join(r1))

	if args.verbose:
		print(f"[weaver] FASTQ written to:\n\t{r1_out}" + (f"\n\t{r2_out}" if args.r2 else ""))
		print(f"[weaver] Total reads: {target_n}")

elif args.command == "extend":
	toolbox.sc_fastq(args.r1, args.r2)

	if args.output:
		os.makedirs(args.output, exist_ok=True)
	
	random.seed(args.seed)
	target_n = args.numReads

	base_r1 = os.path.splitext(os.path.basename(args.r1))[0]
	r1_out = os.path.join(args.output, f"{base_r1}.weaver.fastq")
	toolbox.cp(args.r1, r1_out)

	if args.r2:
		base_r2 = os.path.splitext(os.path.basename(args.r2))[0]
		r2_out = os.path.join(args.output, f"{base_r2}.weaver.fastq")
		toolbox.cp(args.r2, r2_out)

	if args.verbose:
		print(f"[weaver] Copied original reads to:\n\t{r1_out}" + (f"\n\t{r2_out}" if args.r2 else ""))

	with toolbox.smart_open_read(args.r1) as f:
		orig_reads = sum(1 for _ in f) // 4

	if args.verbose:
		print(f"[weaver] Original reads: {orig_reads}")

	if orig_reads >= target_n:
		print(f"[weaver] Requested reads ({target_n}) already satisfied by original")
		sys.exit(0)

	num_to_add = target_n - orig_reads
	if args.verbose:
		print(f"[weaver] Adding {num_to_add} synthetic reads from genome")

	genome_dict = toolbox.read_fasta(args.genome)
	genome = list(genome_dict.values())
	if args.verbose:
		print(f"[weaver] Loaded {len(genome)} sequences from {args.genome}")

	read_len = toolbox.get_read_length(args.r1)
	if args.verbose:
		print(f"[weaver] Read length set to {read_len}")

	with toolbox.smart_open_append(r1_out) as fout1:
		if args.r2:
			with toolbox.smart_open_append(r2_out) as fout2:
				for i in range(num_to_add):
					seq1 = random.choice(genome)

					if len(seq1) < read_len:
						continue

					start = random.randint(0, len(seq1) - read_len)
					syn_seq1 = seq1[start:start + read_len]
					syn_qual1 = "J" * read_len

					seq2 = random.choice(genome)

					if len(seq2) < read_len:
						continue

					start2 = random.randint(0, len(seq2) - read_len)
					syn_seq2 = seq2[start2:start2 + read_len]
					syn_qual2 = "J" * read_len

					hdr = f"WEAVERSYNTH.{i+1} {i+1} length={read_len}"
					fout1.write(f"@{hdr}\n{syn_seq1}\n+{hdr}\n{syn_qual1}\n")
					fout2.write(f"@{hdr}\n{syn_seq2}\n+{hdr}\n{syn_qual2}\n")
		else:
			for i in range(num_to_add):
				seq = random.choice(genome)

				if len(seq) < read_len:
					continue

				start = random.randint(0, len(seq) - read_len)
				syn_seq = seq[start:start + read_len]
				syn_qual = "J" * read_len

				hdr = f"WEAVERSYNTH.{i+1} {i+1} length={read_len}"
				fout1.write(f"@{hdr}\n{syn_seq}\n+{hdr}\n{syn_qual}\n")

	if args.verbose:
		print(f"[weaver] FASTQ written to:\n\t{r1_out}" + (f"\n\t{r2_out}" if args.r2 else ""))
		print(f"[weaver] Reads extended: {num_to_add}")

elif args.command == "synth":
	if args.output:
		os.makedirs(args.output, exist_ok=True)
	
	random.seed(args.seed)
	target_n = args.numReads

	read_len = args.readLength

	genome_dict = toolbox.read_fasta(args.genome)
	genome = list(genome_dict.values())

	if args.verbose:
		print(f"[weaver] Loaded {len(genome)} sequences from {args.genome}")

	r1_out = os.path.join(args.output, f"{args.prefix}_r1.weaver.fastq")
	if args.paired:
		r2_out = os.path.join(args.output, f"{args.prefix}_r2.weaver.fastq")

	if args.verbose:
		print(f"[weaver] Generating {target_n} synthetic read{' pairs' if args.paired else 's'} of length {read_len}")

	with toolbox.smart_open_write(r1_out, False) as fout1:
		if args.paired:
			with toolbox.smart_open_write(r2_out, False) as fout2:
				for i in range(target_n):
					seq1 = random.choice(genome)

					if len(seq1) < read_len:
						continue

					start1 = random.randint(0, len(seq1) - read_len)
					frag1 = seq1[start1:start1 + read_len]
					qual1 = "J" * read_len

					seq2 = random.choice(genome)

					if len(seq2) < read_len:
						continue

					start2 = random.randint(0, len(seq2) - read_len)
					frag2 = seq2[start2:start2 + read_len]
					qual2 = "J" * read_len

					hdr = f"WEAVERSYNTH.{i+1} {i+1} length={read_len}"
					fout1.write(f"@{hdr}\n{frag1}\n+{hdr}\n{qual1}\n")
					fout2.write(f"@{hdr}\n{frag2}\n+{hdr}\n{qual2}\n")
		else:
			for i in range(target_n):
				seq = random.choice(genome)

				if len(seq) < read_len:
					continue

				start = random.randint(0, len(seq) - read_len)
				frag = seq[start:start + read_len]
				qual = "J" * read_len

				hdr = f"WEAVERSYNTH.{i+1} {i+1} length={read_len}"
				fout1.write(f"@{hdr}\n{frag}\n+{hdr}\n{qual}\n")

	if args.verbose:
		print(f"[weaver] Synthetic FASTQ file(s) written to:\n\t{r1_out}" + (f"\n\t{r2_out}" if args.paired else ""))
		print(f"[weaver] Reads synthesized: {target_n}")

else:
	print("\n[weaver] ERROR: subcommand not recognized, choose one of: {reuse, extend, synth}")
	parser.print_help()
