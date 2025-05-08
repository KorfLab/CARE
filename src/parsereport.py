# parsereport

import argparse
import os
import csv
from pathlib import Path


def parse_time_log(filepath):
	"""Parse a single /usr/bin/time -v output file into a dictionary"""
	log = {}
	with open(filepath, "r") as fp:
		for line in fp:
			key, value = line.split(': ')
			log[key.strip()] = value.strip()
	return log


def collect_logs(base_dir):
	"""Walk the results directory and collect parsed logs"""
	results = {}
	for experiment_dir in sorted(Path(base_dir).iterdir()):
		if not experiment_dir.is_dir():
			continue
		experiment = experiment_dir.name
		results[experiment] = {}

		for aligner_dir in sorted(experiment_dir.iterdir()):
			if not aligner_dir.is_dir():
				continue
			aligner = aligner_dir.name
			results[experiment][aligner] = []

			for log_file in sorted(aligner_dir.glob("*.time.log")):
				parsed = parse_time_log(log_file)
				parsed['log_filename'] = log_file.name
				results[experiment][aligner].append(parsed)

	return results


def write_tsv(results, output_dir):
	"""Write out TSV files by experiment"""
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	for experiment in results:
		aligners = results[experiment]
		all_keys = set()
		for aligner in aligners:
			aligner_logs = aligners[aligner]
			for log in aligner_logs:
				for key in log:
					all_keys.add(key)

		fieldnames = sorted(list(all_keys))

		tsv_path = os.path.join(output_dir, experiment + ".tsv")
		f = open(tsv_path, 'w', newline='')
		writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
		writer.writeheader()
		for aligner in aligners:
			logs = aligners[aligner]
			for log in logs:
				writer.writerow(log)
		f.close()

		print("Written", tsv_path)


parser = argparse.ArgumentParser(description="Parse time logs and output TSVs by experiment.")
parser.add_argument('input_dir', 
	help="Path to the results directory")
parser.add_argument('output_dir', 
	help="Path to save TSV files")
args = parser.parse_args()

results = collect_logs(args.input_dir)
write_tsv(results, args.output_dir)
