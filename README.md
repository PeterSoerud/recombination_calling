# Gene conversion and crossing-over events identified from long-read sequencing of primate testes and human sperm

## Overview
This repository contains all the relevant information to reproduce the results and figures in the paper: [Gene conversion and crossing-over events identified from long-read sequencing of primate testes and human sperm](https://www.biorxiv.org/content/10.1101/2024.07.05.601967v1). In this paper, we build a caller that detects recombintion events using long-reads sequencing. Applying a likelihood-based method, we also classify a significant proportion of all reads in each individual, and therefore we can separate the somatic recombinations events from the germline recombination events. Across all individuals in this study, we estimate the gene conversion (GCV) tract length and the amount of GC-biased gene conversion (gBGC).

### Scripts used to estimate the GCV tract lengths and the amount of gBGC:
  [Tract lengths and gBGC](https://github.com/r02ap19/GeneConv)
### Scripts used to classify reads:
  Vinod

### Data from other studies
This repository contains two folders, "tables" and "scritps", where all files in "scripts" are build specificly for this project. However, some files in the "tables" folder are built for other studies, and we list them below:
  - *tables/sd/SD_sort.bed* contains the segmental duplications found in [Vollger et al. 2023](https://www.nature.com/articles/s41586-023-05895-y)
  - *tables/dsb_map/pratto2014.txt* cotains the DSB hotspots found in [Pratto et al. 2014](https://www.science.org/doi/10.1126/science.1256442)
  - *tables/recombination_map/recombination_map_2019.txt* contains the recombination map found in [Halldorsson et al. 2019](https://www.science.org/doi/10.1126/science.aau1043)

## Pipeline
In this section we describe how to run the pipeline that produces the final tables with the classified recombination events.

### Assembly
From raw PacBio HiFi reads, we build a de novo assembly using the software [hifiasm](https://github.com/chhylp123/hifiasm). We run the program with its default parameters, and the code looks like: 
```
hifiasm -o output.asm -t 32 input_reads.fq.gz
```
The program outputs the assembly in a > **_.gfa:_**




