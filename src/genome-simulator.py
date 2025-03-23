import argparse
import random

from toolbox import FTX, anti


def random_exon(n, debug):
	if debug:
		return '[' + '#' * (n-2) + ']'
	return random_seq(n)


def random_flank(n, debug):
	if debug:
		return '~' * n
	return random_seq(n)


def random_intron(n, site, strand, debug):
	if strand == '-':
		site = anti(site)
	if debug:
		return site[0:2] + '-' * (n-4) + site[2:]
	return site[0:2] + random_seq(n-4) + site[2:]


def random_seq(n):
	return ''.join(random.choices('ACGT', k=n))


# CLI
parser = argparse.ArgumentParser(description='Create synthetic fasta & ftx',
								 epilog='structure: ~~~[exon]--intron--[var.exon]--intron--[exon]~~~')
parser.add_argument('name', metavar='<name>',
					help='results in <name>.fa, <name>.ftx')
parser.add_argument('--flank', type=int, default=50, metavar='<int>',
					help='flanking sequence on either side of gene [%(default)i]')
parser.add_argument('--emin', type=int, default=5, metavar='<int>',
					help='minimum middle exon length [%(default)i]')
parser.add_argument('--emax', type=int, default=50, metavar='<int>',
					help='maximum middle exon length [%(default)i]')
parser.add_argument('--estep', type=int, default=1, metavar='<int>',
					help='middle exon step size [%(default)i]')
parser.add_argument('--exon', type=int, default=100, metavar='<int>',
					help='outer exon lengths [%(default)i]')
parser.add_argument('--intron', type=int, default=50, metavar='<int>',
					help='intron lengths [%(default)i]')
parser.add_argument('--chroms', type=int, default=10, metavar='<int>',
					help='number of chromosomes in experiment [%(default)i]')
parser.add_argument('--seed', type=int, metavar='<int>',
					help='set random seed')
parser.add_argument('--double', action='store_true',
					help='create genes on both strands')
parser.add_argument('--noncanonical', action='store_true',
					help='include non-canonical splice sites')
parser.add_argument('--debug', action='store_true')

# NEW ARGS
parser.add_argument('--genomesize', type=int, metavar='<bp>',
					help='Approximate total genome size (bp)')
parser.add_argument('--genesize', type=int, default=500, metavar='<bp>',
					help='Average gene size to estimate genome size [%(default)i]')

arg = parser.parse_args()

# SETUP
if arg.debug:
	arg.flank = 5
	arg.emin = 5
	arg.emax = 6
	arg.exon = 10
	arg.intron = 10
	arg.chroms = 2
	arg.noncanonical = True
	arg.double = True

if arg.seed:
	random.seed(arg.seed)

sites = ['GTAG']
if arg.noncanonical:
	sites.extend(['GCAG', 'ATAC', 'AATT'])
strands = ['+']
if arg.double:
	strands.append('-')

# Determine gene/chromosome counts
if arg.genomesize:
	approx_gene = arg.genomesize // arg.genesize
	arg.chroms = max(1, approx_gene // 50)
	genes_per_chrom = approx_gene // arg.chroms
	print(
		f'Generating ~{approx_gene} genes across {arg.chroms} chromosomes (~{arg.genomesize} bp total)')
else:
	# default behavior
	genes_per_chrom = 10  # adjustable fallback

fa = open(f'{arg.name}.fa', 'w')
fx = open(f'{arg.name}.ftx', 'w')

# MAIN LOOP
for c in range(arg.chroms):
	chrom = f'c{c+1}'
	print(f'>{chrom}', file=fa)
	off = 0
	for g in range(genes_per_chrom):
		for strand in strands:
			for site in sites:
				f1 = random_flank(arg.flank, arg.debug)
				e1 = random_exon(arg.exon, arg.debug)
				i1 = random_intron(arg.intron, site, strand, arg.debug)
				ve = random_exon(random.randint(arg.emin, arg.emax), arg.debug)
				i2 = random_intron(arg.intron, site, strand, arg.debug)
				e2 = random_exon(arg.exon, arg.debug)
				f2 = random_flank(arg.flank, arg.debug)
				e1b = len(f1) + off
				e1e = e1b + len(e1)
				evb = e1e + len(i1)
				eve = evb + len(ve)
				e2b = eve + len(i2)
				e2e = e2b + len(e2)
				exons = ((e1b, e1e-1), (evb, eve-1), (e2b, e2e-1))
				seq = f1 + e1 + i1 + ve + i2 + e2 + f2
				off += len(seq)
				print(seq, file=fa)
				ftx = FTX(chrom, f'g{g}:{site}', strand, exons, '')
				print(ftx, file=fx)

fa.close()
fx.close()
