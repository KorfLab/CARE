# Computational Allocation of Resource Estimate (CARE)

## Todo List

1. Focusing on alignment first
    - Create dev data (obtain small SAM)
    - Develop a genome and read size subsetter by percentage
    - Demo the subsetter with small test data and a tutorial
2. More experiments on alignment performance (momory and time, for multiple aligners)
    - Alignment on various genome size (subsetted genomes)
    - Alignment with various read length (75, 100, 125, 150?)
    - number of reads (1k, 2k, 4k, 8k ...)
    - number of CPUs (1, 2, 4, 8, ...)
3. Determine other tasks to test (qc, trimming, sort, markdup, what else?)
4. Determine programs to test (fastqc, trim-galore, trimmomatic, various alignment tools, samtools and its various commands, ...)
5. Create environment and write push button tester program
6. Test each task with each program in batch
7. Make a figure of resource usage (including time) scaling with task size/complexity (automated)
8. Write a useful manual on cluster computing resource estimation
   - Alignment
   - Other tasks

## Description

This project aims to test resource allocation and its effects on the performance of various bioinformatics tools on computing clusters. This effort should enable ultra fast resource estimate for most future tasks performed on a cluster, whether it is managed through Slurm or not.

Multiple figures of each tool shall be made to provide visual showcase of performance (cpu time) versus resource provided.

A user manual should be written at the end of the project to teach future users how cluster computing works and how to correctly and accurately identify the right amount of resources needed for a task given the size of the task itself.

An html based UI may be created for fast estimation, but that is the least priority comparing to everything else.
