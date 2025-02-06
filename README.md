# Gene conversion and crossing-over events identified from long-read sequencing of primate testes and human sperm

## Overview
This repository contains all the relevant information to reproduce the results and figures in the paper: [Gene conversion and crossing-over events identified from long-read sequencing of primate testes and human sperm](https://www.biorxiv.org/content/10.1101/2024.07.05.601967v1). In this paper, we build a caller that detects recombintion events using long-reads sequencing. Applying a likelihood-based method, we also classify a significant proportion of all reads in each individual, and therefore we can separate the somatic recombinations events from the germline recombination events. Across all individuals in this study, we estimate the gene conversion (GCV) tract length and the amount of GC-biased gene conversion (gBGC).

### Scripts used to estimate the GCV tract lengths and the amount of gBGC:
  [Tract lengths and gBGC](https://github.com/r02ap19/GeneConv)
### Scripts used to classify reads:
  [Classification of PacBio HiFi reads](https://github.com/vinodsinghjnu/Proj_GermReadsClassification/blob/main/README.md)

### Data from other studies
This repository contains two folders, "tables" and "scritps", where all files in "scripts" are build specificly for this project. However, some files in the "tables" folder are built for other studies, and we list them below:
  - `tables/sd/SD_sort.bed` contains the segmental duplications found in [Vollger et al. 2023](https://www.nature.com/articles/s41586-023-05895-y)
  - `tables/dsb_map/pratto2014.txt` cotains the DSB hotspots found in [Pratto et al. 2014](https://www.science.org/doi/10.1126/science.1256442)
  - `tables/recombination_map/recombination_map_2019.txt` contains the recombination map found in [Halldorsson et al. 2019](https://www.science.org/doi/10.1126/science.aau1043)

## Pipeline
In this section we describe how to run the pipeline that produces the final tables with the classified recombination events.

### Assembly
From raw PacBio HiFi reads, we build a de novo assembly using the software [hifiasm](https://github.com/chhylp123/hifiasm). We run the program with its default parameters, and the code looks like: 
```
hifiasm -o output.asm -t 32 input_reads.fq.gz
```
The program outputs the assembly in a `.gfa` file, and to convert the file into fasta format, we run the code below. To complete the assembly we also build an index file using *samtools*.
```
awk '/^S/{{print ">" $2;print $3}}' gfa_file_in | fold > fasta_out
samtools faidx fasta_in > index_out
```
### Mapping
One of the most crucial steps in this study is mapping the reads back against the de novo assembly. For this task we use the software [pbmm2](https://github.com/PacificBiosciences/pbmm2) with the settings specified below:
```
pbmm2 index fa_in index_out
pbmm2 align fa_in fastq_in bam_out --unmapped --preset HIFI --sort -j 12 -J 4 -m 8G -c 99 -l 2740
```
Since we by default would reject all events found on soft-clipped reads, we remove these before further analysis. The soft-clippede reads are removed using *samtools*:
```
samtools view -h bam_in | awk '$6 !~ /H|S/{{print}}' | samtools view -bS > bam_no_clips
samtools index bam_no_clips
```
### Recombination calling
Before calling the recombination events, we run `scripts/calling/coverage.py` to obtain the coverage distribution of all contigs great than 1 Mb. The script take a bam file as input and outputs a table with counts of the coverage. This table is used as input to `scripts/calling/coverage.Rmd`, which calculates the upper and lower limit used as threshold for the recombination callling. 

We then run `scripts/calling/recombination.py`, which imports functions from `scripts/calling/utils.py`. The recombination caller takes a fasta, fasta.fai, and the bam file as input files. It also takes several thresholds as input, and these numbers can easily be changed in the script. The script outputs six different files:
  - A tab-separated file containing all the called recombination events
  - A tab-separated file containing all the recombination calls that were discarded
  - A tab-separated file containing the inter-SNV distances
  - A tab-separated file containing all the SNV calls that were discarded
  - A tab-separated file containing the positions of runs of homozygosity (>1e5 bps between SNVs)
  - A fastq file containing alle the reads in which an event was called

### Reference genomes
The reference genomes used in this study are listed below:
  - [Homo sapiens](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/)
  - [Pan troglodytes](https://genomeark.s3.amazonaws.com/index.html?prefix=species/Pan_troglodytes/mPanTro3/assembly_curated/)
  - [Gorilla gorilla](https://genomeark.s3.amazonaws.com/index.html?prefix=species/Gorilla_gorilla/mGorGor1/assembly_curated/)
  - [Hylobates pileatus](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_021498465.1/) (Sequenced gibbon is a Hylobates lar)
  - [Papio papio](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_028645565.1/)
  - [Macaca nemestrina](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_043159975.1/)
    
To convert the contig positions of the called events to positions in the corresponding reference genome, we map the reads containing recombination events using *pbmm2*. We then use *samtools* to extract the relavant information from the mapping. The code is shown below:
```
pbmm2 align reference_genome.fa fastq_in bam_out --preset HIFI --sort -j 12 -J 4 -m 8G
samtools view bam_in | awk -v OFS='\t' '{{print $1,$2,$3,$4}}' > bed_out
```
Since we also want to know all events in human coordinates, we also map all non-human samples against the human reference genome (T2T-CHM13v2.0)

### Classify reads with events
We are capable of classifying a significant amount of the reads from each individual, and to add the classification to the reads with a called recombination event, we simply intersect the two data tables. However, to join in bash we need to sort the tables by read name:
```
awk 'NR==1{print; next} {print | "sort -k5,5"}' recombination_table.bed > recombination_table_sorted.bed
awk 'NR==1{print; next} {print | "sort -k1,1"}' classification_table.tsv > classification_table_sorted.tsv
join -1 5 -2 1 --header recombination_table_sorted.bed classification_table_sorted.tsv > recombination_plus_methylation.tsv
```

### Process recombination events
At this stage of the analysis, we have a set of candidate recombination events for each individual. However, we cannot classify the reads for the Guinea baboon, Lar Gibbon, and pig-tailed macaque because of their evolutionary distance to humans. Therefore we run `/scripts/process_calls/split_raw_calls_unclassified.Rmd` for these three species, and we run `/scripts/process_calls/split_raw_calls_classified_sperm.Rmd` for the rest of the samples. The script outputs four tables, one for each of the recombination types introduced in the paper. 

### Manual curation
For each sample we manual curate each of its four recombination classes. The results of the curation can be found in `/tables/manual_curation`.

### Final files
Lastly, we run `/scripts/process_calls/split_raw_calls_classified_sperm.Rmd` and `/scripts/process_calls/split_raw_calls_classified_sperm.Rmd`, which uses the files from the manual curation to filter the recombination events. The output files from these scripts constitute the basis of all further analyses. 



