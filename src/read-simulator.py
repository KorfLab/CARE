import argparse
import random
import sys

from toolbox import FTX, generator, anti


def generate_reads(gftx, chrom, size):
    # create indexes
    dna = []  # dna positional index
    rna = []  # rna sequence
    for beg, end in gftx.exons:
        for i in range(end - beg + 1):
            coor = i + beg
            dna.append(coor)
            rna.append(chrom[coor])

    # generate reads and their ftx annotations
    for i in range(len(rna) - size + 1):
        coor = [dna[i + j] for j in range(size)]
        exons = []
        beg = coor[0]
        seen = 0
        for j in range(size - 1):
            d = coor[j + 1] - coor[j]
            if d > 1:
                end = beg + j - seen
                exons.append((beg, end))
                seen += end - beg + 1
                beg = coor[j + 1]
        exons.append((beg, beg + j - seen + 1))
        rftx = FTX(gftx.chrom, gftx.name, gftx.strand, exons, 'r')
        read = ''.join([chrom[beg:end + 1] for beg, end in exons])

        yield rftx, read


parser = argparse.ArgumentParser()
parser.add_argument('fasta', help='fasta file')
parser.add_argument('ftx', help='ftx file')
parser.add_argument('--readlength', type=int, default=100, metavar='<int>',
                    help='[%(default)i]')
parser.add_argument('--numreads', type=int, metavar='<int>',
                    help='Desired total number of reads (approx)')
parser.add_argument('--seed', type=int, default=0, metavar='<int>',
                    help='set random seed')
parser.add_argument('--double', action='store_true',
                    help='produce reads from both strands')

arg = parser.parse_args()

if arg.seed != 0:
    random.seed(arg.seed)

genes = 0
reads = 0
bases = 0
all_reads = []

# First pass: collect possible reads
for cname, cseq, gtfxs in generator(arg.fasta, arg.ftx):
    for gftx in gtfxs:
        for rftx, rseq in generate_reads(gftx, cseq, arg.readlength):
            all_reads.append((rftx, rseq, '+'))
            if arg.double:
                all_reads.append((rftx, anti(rseq), '-'))

total_possible = len(all_reads)
print(f'Total possible reads: {total_possible}', file=sys.stderr)

# Determine downsampling
if arg.numreads and arg.numreads < total_possible:
    sampled_reads = random.sample(all_reads, arg.numreads)
else:
    sampled_reads = all_reads

# Output reads
for rftx, rseq, strand in sampled_reads:
    print(f'>{rftx}{strand}')
    print(rseq)
    reads += 1
    bases += len(rseq)

print(f'genes: {genes}', f'reads: {reads}', f'bases: {bases}',
      sep='\n', file=sys.stderr)
