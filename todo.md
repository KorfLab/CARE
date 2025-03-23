# Computational Allocation of Resource Estimate (CARE)

## Todo List

1. Make test data
    1a. Re-learn Slurm commands and job submission processes
    1b. Determine tasks to test (qc, trimming, read alignment, sort, markdup, what else?)
2. Determine programs to test (fastqc, trim-galore, trimmomatic, various alignment tools, samtools and its various commands, what else?)
3. Create environment and write push button tester program
4. Test each task with each program in batch
5. Make a figure of resource usage (including time) scaling with task size/complexity
6. Write a useful manual on cluster computing resource estimation
7. Make a push button UI for quick estimate

## Description

This project aims to test resource allocation and its effects on the performance of various bioinformatics tools on computing clusters. This effort should enable ultra fast resource estimate for most future tasks performed on a cluster, whether it is managed through Slurm or not.

Multiple figures of each tool shall be made to provide visual showcase of performance (cpu time) versus resource provided.

A user manual should be written at the end of the project to teach future users how cluster computing works and how to correctly and accurately identify the right amount of resources needed for a task given the size of the task itself.

An html based UI may be created for fast estimation, but that is the least priority comparing to everything else.
