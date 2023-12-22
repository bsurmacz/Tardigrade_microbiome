# Tardigrade_microbiome
 
 

#### This repository contains the complete pipelines for the analysis and visualization of 16S rRNA V4 amplicon sequencing data for the analysis of tardigrade microbiome and re-analysis of data from other tardigrade microbiome studies.  Libraries were sequenced across three 2 x 300 bp Illumina MiSeq lanes, in multiplex with samples from other projects


### The starting data:

#### Experiment 1: Broad preliminary survey
#### Experiment 2: Comparing replicate cultures
#### Experiment 3: Comparing sample preparation approaches 
#### Experiment 4: Species comparisons

The workflow comprised 4 distinct steps. They are listed below, and either discussed in more details, or linked to further down.

1. Assembly and extraction of project-specific libraries from multiplexed sequencing datasets and merging into one dataset 
2. Basic analysis of the amplicon data: filtering, denoising, OTU picking
3. Decontamination
4. Splittinng the data from separate experiments - to maintain OTU IDs for the whole study,

## 1. Extraction of project-specific libraries from multiplexed sequencing datasets and merging into one dataset 

## 2. Basic analysis of the amplicon data: filtering, denoising
For this step we need the following: 
- [vsearch v2.15.2_linux_x86_64](https://github.com/torognes/vsearch)
- [usearch v11.0.667_i86linux32](https://www.drive5.com/usearch/)


## 3. Splittinng the data from separate experiments
As we wanted to run the analyses (decontamination, filtering) separately for different experiments (maintaining the OTUs' names across experiments), we need to separate the data originating from different experiment.

## 4. Decontamination of the data in each experiment separately
For this step we need the following: 
- Python 3.7.8
- custom decontamination script **decontaminate.py**


The outputs of this script are:
- **Table_with_classes.txt** - where every zOTU is assigned to Symbiont, Other, PCR or Extraction Contaminant and PCR or Extraction Spikein class
- **Statistics_table.txt** - with statistics about every library composition in terms of e.g contamination or spikein percentage 
- **Decontaminated_zOTU_table.txt** - where all contaminants and spikeins are deleted, as well as libraries that sum of those were higher than ThresholdC
- **Decontaminated_OTU_table.txt** - table based on Decontaminated zOTU table and otus.tax file








