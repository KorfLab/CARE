# Subsetter Tutorial: Shrinking a Genome and Reads to Specified Percentage

This tutorial demonstrates how to use the `subset.py` script under the `/src` directory to shrink both a genome and a set of reads based on a SAM file. This can be useful for generating scaled-down test datasets for performance benchmarking or debugging while ensuring the scaled down reads would all align to the scaled down genome.

We will walk through shrinking both a genome FASTA and corresponding FASTQ reads to 50% of their original capacity. 

Here a very small dataset is used:

- A partial *C. elegans* genome: `ref/1pCe.fa` (1% of the original)
- A set of reads known to align within this subset: `readsx3000.fastq`

For showcasing purposes, we will be using `STAR` as the sequence aligner to obtain the SAM file. Please make sure you have `STAR` installed in your working environment.

---

## Step 1: Index the Genome

First, index the genome with STAR. Since this is a small genome (1%), we set `--genomeSAindexNbases` to 8.

```bash
STAR --runMode genomeGenerate \
     --runThreadN 4 \
     --genomeDir ref/ \
     --genomeFastaFiles ref/1pCe.fa \
     --genomeSAindexNbases 8
```

---

## Step 2: Align Reads to the Genome

Use STAR in `alignReads` mode to align the reads and generate a SAM file.

```bash
STAR --runMode alignReads \
     --runThreadN 4 \
     --genomeDir ref/ \
     --readFilesIn readsx3000.fastq \
     --outFileNamePrefix star_output/ \
     --outSAMunmapped Within
```

This will create a SAM file at `star_output/Aligned.out.sam`.

---

## Step 3: Shrink Genome and Reads Using `subset.py`

Now we use `subset.py` to shrink the genome and reads based on the alignment information. With `-p 0.5` or `--scale 0.5`, it will create a new genome containing only the first 50% of each chromosome, and **only retain reads that still fully align** within this shortened genome.

Make sure you are in the /demo directory.

```bash
python3 ../src/subset.py \
    --genome ref/1pCe.fa \
    --sam star_output/Aligned.out.sam \
    --scale 0.5 \
    --r1 readsx3000.fastq \
    -v
```

This will produce:

- `1pCe.shrunk.fa`: a FASTA file with 50% of each original sequence, within the /ref directory
- `readsx3000.shrunk.fastq`: a FASTQ file with only the reads that still align to the shrunk genome

You will also get printed output summarizing how many reads were kept and the new sizes of the genome and read files.

If you prefer gzipped output, add the `-z` flag will create the following file instead:

- `readsx3000.shrunk.fastq.gz`: compressed FASTQ file with only the reads that still align to the shrunk genome

Note `STAR` won't take a compressed FASTQ file as input during read alignment. Use `gzip -dk readsx3000.shrunk.fastq.gz` to uncompress it first.

That's it. You now have a downsized genome and corresponding reads, which can be used for testing or benchmarking on smaller scales.

---

## Step 4: Benchmark STAR Performance with and without Unaligned Reads

One goal of this subsetter is to create a smaller genome with a corresponding set of filtered reads such that **every read aligns** to the shrunk genome. This presents an ideal alignment scenario for an aligner.

To understand how **unaligned reads impact aligner performance**, we will perform a second alignment using the **original (unfiltered)** reads against the **shrunk genome**, where many reads are expected to be unaligned. To keep the setup uniform, we force the number of threads to 4 for both runs.

We will use `/usr/bin/time` to measure alignment time in both scenarios.

---

### 4.1 Align Shrunk (Filtered) Reads to Shrunk Genome

This is the "ideal" case: every read aligns successfully.

Make sure to index the shrunk genome first.

```bash
mkdir ref/shrunk
STAR --runMode genomeGenerate \
     --runThreadN 4 \
     --genomeDir ref/shrunk \
     --genomeFastaFiles ref/1pCe.shrunk.fa \
     --genomeSAindexNbases 8
```

```bash
/usr/bin/time -v STAR --runMode alignReads \
    --runThreadN 4 \
    --genomeDir ref/shrunk \
    --readFilesIn readsx3000.shrunk.fastq \
    --outFileNamePrefix star_shrunk_to_shrunk/ \
    --outSAMunmapped Within
```

You’ll see timing statistics such as:
- Elapsed wall clock time
- User and system CPU time
- Maximum resident set size (memory)

Take note of these numebrs as we will compare the results to those from the next scenario.

---

### 4.2 Align Original (Unfiltered) Reads to Shrunk Genome

In this case, many reads will fail to align, which will increase runtime.

```bash
/usr/bin/time -v STAR --runMode alignReads \
    --runThreadN 4 \
    --genomeDir ref/shrunk \
    --readFilesIn readsx3000.fastq \
    --outFileNamePrefix star_unfiltered_to_shrunk/ \
    --outSAMunmapped Within
```

This simulates a scenario where half of the reads won't map to the target genome.

---

## Summary

| Alignment Run                  | Expected Alignment Rate | Purpose                             | Runtime |
| ------------------------------ | ----------------------- | ----------------------------------- | ------- |
| Filtered reads - Shrunk genome | =100%                   | Ideal-case performance benchmark    | ????? s |
| Original reads - Shrunk genome | <100%                   | Stress-test for handling mismatches | ????? s |

Use the timing results to determine whether handling unaligned reads adds a measurable performance cost to the aligner.

Hint: the difference is **HUGE**.
