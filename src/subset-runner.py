import argparse
import os
import sys
import yaml


#############
# functions #
#############

def run_cmd(cmd):
	"""Run a command and exit when failed"""
	print(f"[run cmd] {cmd}")
	if os.system(cmd) != 0:
		sys.exit(f'FAILED: {cmd}')


############
# argparse #
############

parser = argparse.ArgumentParser(
	description="Batch subset genome using percentages from YAML")
parser.add_argument("-g", "--genome", required=True, 
	help="Genome FASTA file")
parser.add_argument("-y", "--yaml", required=True,
	help="YAML file of a list of percentages")
parser.add_argument("-a", "--aligner", nargs='+',
	help="Aligners to prepare genome folders for")

args = parser.parse_args()


########
# main #
########

gabs = os.path.abspath(args.genome)
gdir = os.path.dirname(gabs)
gbsn = os.path.splitext(os.path.basename(gabs))[0]

with open(args.yaml) as ymlin:
	pcts = yaml.safe_load(ymlin)

for p in pcts:
	if not isinstance(p, float) or not (0 < p <= 1):
		print(f"Invalid percentage in YAML: {p}, must be a float in (0, 1]")
		continue

	if p == 1:
		dst = os.path.join(gdir, "ref-100p.fa")
		cmd = ["cp", "-f", gabs, dst]
		run_cmd(' '.join(cmd))
		continue

	cmd = [
		"python3", "subset.py", "xfa",
		"-g", args.genome,
		"-p", f"{p:.2f}",
		"-v"
	]
	run_cmd(' '.join(cmd))

	src = os.path.join(gdir, f"{gbsn}.shrunk.fa")
	dst = os.path.join(gdir, f"ref-{int(p * 100)}p.fa")
	cmd = ["mv", "-f", src, dst]
	run_cmd(' '.join(cmd))


########
# Move #
########

if args.aligner:
	print(f"\n[INFO] Organizing genome subsets for aligners: {args.aligner}")

	refs = [f for f in os.listdir(gdir) if f.startswith("ref-") and f.endswith("p.fa")]

	for aligner in args.aligner:
		aligner_dir = os.path.join(gdir, f"genomes-{aligner}")
		os.makedirs(aligner_dir, exist_ok=True)

		for ref in refs:
			src = os.path.join(gdir, ref)
			dst = os.path.join(aligner_dir, ref)
			cmd = ["cp", "-f", src, dst]
			run_cmd(' '.join(cmd))
