# Design Outline

Overall container - Conda

## Processes

- QC of reads
- Trimming of low quality reads and markers
- Read alignment
- Read sorting
- Read duplication marking
- Go further until variant calling?
- Anything missing?

## Programs

### QC

- FastQC
- MultiQC
- 

### Trimming

- Trimmomatic
- Trim-galore
- 

### Alignment

- bwa mem
- star
- hisat2
- minimap2
- 

### Samtools Various processes

- index
- sort
- markdup
- fixmate
- 

### Designs

- goal: optimal CPU time with minimal resoures of various task sizes
- method: repetative tests allow graph of time versus resources
- task size: 
  - File sizes
    - number of reads
    - read length
    - genome size
- Resources:
  - Memory allowance
  - CPU available
  - anything else?

### End goal

Graphs of task CPU time versus each resource

Hypothesis: Should see a peak and diminishing return with 
increased resources for tasks

A tool that can estimate optimal resource allocation given
a job size (number of reads, genome size, read length, 
number of tasks)
