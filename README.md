This repository contains all the software and scripts used for the following study.

Garcia-Nieto PE, Morrison AJ, Fraser HB. The somatic mutation landscape of the human body. BioRxiv. 2019 [link](https://www.biorxiv.org/content/10.1101/668624v1)
Garcia-Nieto PE, Morrison AJ, Fraser HB. The somatic mutation landscape of the human body. Genome Biology. 2019 [link](https://www.biorxiv.org/content/10.1101/668624v1)

The contents of this repository encompass:
- A method to call DNA somatic mutations from the GTEx v7 RNA-seq data:
    - Mapping raw reads to the human genome Hg19.
    - Calling mutations.
    - Identifying and removing false-positive mutations
- All analysis included in the study on the mutations:
    - Method validation using exome sequencing data
    - General characterization of the landscape of mutations.
    - Phenotypic associations.
    - Molecullar associations -- chromatin and gene expression.
    - Strand asymmetry analysis.
    - Mutation
    - Selection dynamics.
    - Assessment of cancer-like characteristic of somatic mutations
    
For a complete description of the methods and the results of theses analysis please refer to the publication

## Requirements
- Unix -- all software was run on bash
- snakemake -- for all pipelines
- BCFtools 1.6
- bedtools 2.26
- samtools 1.6
- R 3.4
    - Tidiverse
    - Bioconductor
        - GenomicRanges
- python 3.6
    - pysam
    - pandas
    - numpy

## Overview

All raw RNA reads were downloadded from dbGaP using GTEx v7 data. 

The execution outline is as follows:

1. Setting up configuration file and downloading all required auxiliary files.
2. Mapping to the human genome.
3. Creation of readp pileup files for positions contating two different sequence calls.
3. Per-position addition of potential artifacts for the alternate allele.
4. Mutation calling and removal of likely false-positives.
5. Elimination of potential systematic artifacts (Panel of Normals).
6. Elimination of hyper-mutated samples.
7. High-level analyses

An in-depth description on how to execute each of these steps is described below.

## Seting up configuration file and download auxiliary data

### Downloading auxiliary files

Follow this [link](https://drive.google.com/a/stanford.edu/file/d/1v9ZIfkMmi7q8yh_lkn2BHFB4XDYxsERx/view?usp=sharing)

Uncompress, and move to a static location. These files are a variety of genome file, gene annotations, 

### Configuring

Modify the file `config.json`, all buckets should be self explanatory. The file is setup as it was used for the original study, which was run on the Sherlock cluster at Stanford University

It's worth noting the following:

- "genomeDir" is the path for the genome STAR index (see mapping for details in how to make it)
- Any bucket containing '\*/auxiliaryFiles/[...]' should be replace with the location of the downloaded and uncompress file above: 'new\_path/auxiliaryFiles/[...]'
- "vcfFile" and "vcfFileCommon" are paths to the genotype vcf files from the GTEx project, please access those through dbGaP, feel free to contact me for any guidance on this.
- "projectDir" and "scratchDir" should be identical, this the root path where all results will be stored.
- Most other buckets can be left unmodified

#### Cluster and parallelization

Most pipelines are designed to run many jobs in parallel. Included there are some extra files and scripts to run snakemake on a SLURM cluster:


- The file 'cluster.json' contains rule-specific specifications and can be used on snakemake.
- submit.py is a custom submission script to be used for job submission on snakemake.
- jobState is a bash script to be used by snakemake to check the sate of jobs, **this script should be added to the bash `$PATH`.**

