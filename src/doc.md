# subset Documentation

This program trims a genome FASTA to a given percentage and/or filters FASTQ reads based on whether they still align after genome truncation. Useful for downscaling datasets for benchmarking, debugging, or running lightweight tests with full mapping logic intact.

There are two modes, specified via subcommands:

- `xfa`: shrink the genome only
- `xfq`: shrink the genome and filter reads based on a SAM file

---

## Basic Usage

```bash
python3 subset.py <subcommand> [options]
```

Subcommands:

- `xfa` — cut genome only
- `xfq` — cut genome and filter reads based on alignment

---

## Subcommand: `xfa`

Trim each chromosome in a genome FASTA to the first X% of its length.

### Required

- `-g`, `--genome`: specify input FASTA file
- `-p`, `--scale`: float in (0, 1], percent of each sequence to keep

### Optional

- `-v`, `--verbose`: print progress and file size info

### Behavior

- Writes a new FASTA file named `[input].shrunk.fa`.
- If `-p`, `--scale` is `1.0` or `1` when using this subcommand, no file is written, and the script exits immediately.

---

## Subcommand: `xfq`

Given a SAM file and FASTQ reads, keep only the reads that still align within the shrunken genome.

### Required

- `-g`, `--genome`: specify input FASTA file, original or shrunk
- `-s`, `--sam`: specify alignment file in SAM format
- `-p`, `--scale`: float in (0, 1]; if 1.0, genome is only read and left unchanged
- `--r1`: FASTQ file (single-end or first pair)

### Optional

- `--r2`: second FASTQ file for paired-end mode
- `-z`, `--gzip`: gzip the output FASTQ files
- `-v`, `--verbose`: print filtering stats and file size info

## Behavior

- Writes new FASTQ files named `[r1].shrunk.fastq` and optionally `[r2].shrunk.fastq`.
- If `-z` is set, output files will have the `.gz` extension.
- If `-p`, `--scale` < 1.0:
  - Genome is shrunk and saved as `[input].shrunk.fa`.
  - Reads are filtered based on alignment to the shrunk genome.
- If `-p`, `--scale` == 1.0:
  - Genome is not modified and no new FASTA file is created.
  - Reads are filtered using original alignment to the full genome.
- Reads that fully align within the genome boundary will be kept.
- For paired-end reads, both mates must fully align to be kept.

---

## Examples

Shrink genome to 25%:

```bash
python3 subset.py xfa -g genome.fa -p 0.25
```

Shrink genome to 50% and keep only reads still aligning; show extra progress info:

```bash
python3 subset.py xfq -g genome.fa -s aligned.sam -p 0.5 --r1 reads.fq -v
```

Only filter mapped paired-end reads without shrinking the genome:

```bash
python3 subset.py xfq -g genome.fa -s aligned.sam -p 1.0 --r1 r1.fq --r2 r2.fq
```

---

## Notes

This program:
- always writes outputs to the same directory as the inputs
- filters reads based on CIGAR and end position vs. genome sequence end position
- assumes input SAM uses the same genome as specified in `--genome`
