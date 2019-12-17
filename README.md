This repository contains all the software and scripts used for the following study.

**Garcia-Nieto PE, Morrison AJ, Fraser HB. The somatic mutation landscape of the human body. BioRxiv. 2019. [link](https://www.biorxiv.org/content/10.1101/668624v1)**

**Garcia-Nieto PE, Morrison AJ, Fraser HB. The somatic mutation landscape of the human body. Genome Biology. 2019. [link](NA)**

The contents of this repository encompass:
- A method to call DNA somatic mutations from the GTEx v7 RNA-seq data:
    - Mapping raw reads to the human genome Hg19.
    - Calling mutations.
    - Identifying and removing false-positive mutations
- All analyses included in the study:
    - Method validation using exome sequencing data
    - General characterization of the landscape of mutations.
    - Phenotypic associations.
    - Molecular associations -- chromatin and gene expression.
    - Strand asymmetry analysis.
    - Selection dynamics.
    - Assessment of cancer-like characteristics of somatic mutations
    
***This document is not intended to be a thorough description of the methods or results. It is a guide that serves as a reproducibility reference for the code use in the aforementioned study.***

For a complete description of the methods and results please refer to the publication.

## Requirements
- Unix -- all software was run on bash in Ubuntu 14.04
- snakemake -- for all pipelines
- STAR 2.5.2a -- sequence aligner
- kent from ucscGenomeBrowser
- BCFtools 1.6
- bedtools 2.28
- samtools 1.6
- R 3.4:
    - MASS
    - MatrixEQTL
    - Only output alignments in library
    - Rsamtools
    - Rtsne
    - biomaRt
    - cluster
    - dndscv
    - ggplot2
    - ggrepel
    - gplots
    - gtools
    - metap
    - proto
    - rSubmitter
    - reshape
    - rjson
    - seqinr
    - topGO
    - Tidyverse
    - Bioconductor:
        - GenomicRanges
        - Biostrings
        - GenomicFeatures
        - GenomicRanges
        - Homo.sapiens
        - SNPlocs.Hsapiens.dbSNP144.GRCh37
        - org.Hs.eg.db
- python 3.6
    - pysam
    - pandas
    - numpy
    - scipy

## Overview

All raw RNA reads were downloaded from dbGaP using GTEx v7 data. 

The execution outline is as follows:

1. Setting up configuration files and downloading all required auxiliary files.
2. Creating a genome index for the STAR sequence aligner.
3. Mapping to the human genome.
4. Creation of read pileup files for positions having two different sequence calls.
5. Per-position addition of potential artifacts for the alternate allele.
6. Mutation calling and removal of likely false-positives.
7. Elimination of potential systematic artifacts (Panel of Normals).
8. Elimination of hyper-mutated samples.
9. High-level analyses.

An in-depth description on how to execute each of these steps is described below.

## Setting up 

#### Downloading auxiliary files

Follow this [link](https://drive.google.com/a/stanford.edu/file/d/1v9ZIfkMmi7q8yh_lkn2BHFB4XDYxsERx/view?usp=sharing)

Uncompress it, and move it to a static location. These files are a variety of genome files, gene annotations, data from external databases (e.g. COSMIC, Roadmap epigenomics, TCGA), RNA edit info, etc.

#### Configuring

Modify the file `config.json`, all buckets should be self-explanatory. The file is setup as it was used for the original study, which was run on the Sherlock cluster at Stanford University.

It's worth noting the following:

- `genomeDir` is the path for the genome STAR index (see mapping for details in how to make it)
- Any bucket containing `\*/auxiliaryFiles/[...]` should be replaced with the location of the downloaded and uncompressed file [above](#downloading-auxiliary-files): `new\_path/auxiliaryFiles/[...]`
- `vcfFile` and `vcfFileCommon` are paths to the genotype vcf files from the GTEx project, please access those through dbGaP, feel free to contact me for any guidance on this.
- `projectDir` and `scratchDir` should be identical, this the root path where all results will be stored.
- Most other buckets can be left unmodified

#### Cluster and parallelization set up

Most pipelines are designed to run manextey jobs in parallel. Included there are some extra files and scripts to run snakemake on a SLURM cluster:


- The file 'cluster.json' contains rule-specific specifications and can be used on snakemake.
- submit.py is a custom submission script to be used for job submission on snakemake.
- jobState is a bash script to be used by snakemake to check the sate of jobs, **this script should be added to the bash `$PATH`**

## Creating genome index for STAR

This is the only step not included as a download or a step in any pipeline due to large memory and disk usage.

#### Inputs
- Human genome Hg19 with common SNPs masked as "N". Included in the `auxiliaryFiles`, the path to the file should be the `genomeFasta` bucket in `config.json`
- Genome annotation in GTF format. Included in the `auxiliaryFiles`, the path to the file should be the `auxiliaryFiles:GTF_gtex` bucket in `config.json`

#### Outputs
- Genome index for STAR

#### Execution
Replace variables ($) with the associated value in the `config.json` file

```bash
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ${genomeDir} --genomeFastaFile ${genomeFasta} --sjdbGTFfile ${auxiliaryFiles:GTF_gtex} --sjdbOverhang 99
```

## Mapping raw RNA-seq reads and creating pileup files

The following pipeline includes:
- Mapping to the human genome, sorting and indexing bam files, removing duplicate reads.
- Creating coverage maps.
- Creating read pileup files for positions covered by two different sequence calls.
- Adding per-site bias information to the pileup files: read position bias, mapping quality bias, sequence quality bias, strand bias, variant distance bias. 

A full description of these biases can be found in the Methods section of the publication.

#### Inputs
- Raw rna-seq paired fastq files. These should be downloaded from dbGaP GTEx v7, please contact me for support on this. The path to these files should be the `fastqDir` bucket in `config.json`. File names should have the format `{sample}_1.fastq.gz` `{sample}_2.fastq.gz`
``nformation. Included in the `auxiliaryFiles` download, the path to the file should be the `sraTable` bucket in `config.json`
- Genome size table. Included in the `auxiliaryFiles` download, the path to the file should be the `genomeSizeFile` bucket in `config.json`
- A variety supporting files included in `auxiliaryFiles`

#### Outputs
- Sorted, bam files without duplicate reads. Files will be located in `$mappingDir/{tissue}/{sample}_RmdupSortedAligned.out.bam`, where `$mappingDir` is the `mappingDir` bucket in `config.json`
- Coverage maps `$root/depth_bam/{sample}.bed.gzip`. Where `$root` is the `projectDir` bucket in `config.json`
- Pileup files for positions with two sequence calls `$pileupDir/{tissue}/{sample}.txt`. Where `$pileupDir` is relative and is the `pileupDir` bucket in `config.json`

#### Execution

Linear execution (this is not feasible as it would take a very large time to finish):
```bash
cd mappingMutationCalling_pipeline
``oing --restart-times 2 
```

Parallel execution in a SLURM cluster:
```bash
cd mappingMutationCalling_pipeline
snakemake --snakefile SnakefilePileup.smk --restart-times 2 --keep-going --max-jobs-per-second 3 --max-status-checks-per-second 0.016 --cluster-config ../cluster.json --cluster-status jobState --jobs 500 --keep-going --cluster "../submit.py"
```

## Mutation calling

The following pipeline includes:
- Calling mutations based on coverage and number of reads supporting alt allele, while incorporating sequencing error probabilities.
- Labelling and removing false positives:
    - Black-listed regions.
    - RNA edits.
    - Splicing junction errors.
    - Sequencing errors.
    - All biases described in the previous section.
- Remove mutations from systematic artifacts based on a pseudo Panel of Normals (PoN).
- Remove hypermutated samples.

A full description of the false-positive removal process can be found in the Methods section of the publication.

#### Inputs
``o sequence calls `$pileupDir/{tissue}/{sample}.txt`. *From the previous section.*
- A variety of supporting files included in `auxiliaryFiles`

#### Outputs
- Mutation maps per sample `${mutationCountDir:map}/{tissue}/n6_0.0_0.7/{sample}.txt`.
- Counts of mutations following a 6-type profile (C\>T, C\>G, C\>T, C\>A, T\>G, T\>A, T\>C) `$mutationCountDir:$context/{tissue}/n6_0.0_0.7/{sample}.txt`
- Counts of mutations with pentanucleotide context per sample `$mutationCountDir:$context/{tissue}/n6_0.0_0.7/{sample}.txt`

#### Execution

Linear execution:
```bash
cd mutationCalling
snakemake --snakefile SnakefilePileup.smk --printshellcmds --keep-going --restart-times 2 
```

Parallel execution in a SLURM cluster:
```bash
cd mutationCalling
snakemake --restart-times 2 --keep-going --max-jobs-per-second 3 --max-status-checks-per-second 0.016 --cluster-config ../cluster.json --cluster-status jobState --jobs 500 --keep-going --cluster "../submit.py"
```

## General analyses

``lysis and plotting scripts for:
- Correction of mutations based on sequencing depth.
- Artifact analysis.
- Mutation statistics -- variant allele frequency, mutation density across subjects and samples, etc.
- Cross-tissue comparison of mutation maps.
- Phenotypic associations (age, ethnicity, sex).
- Mutation profile analysis (tSNE).
- Gene expression association with mutation load.
- Association between gene expression and mutation load.

These analyses span Figures 1,2,4 in the publication.

##### Inputs
- Mutation maps and counts from [above](#mutation-calling)
- A variety of supporting files included in `auxiliaryFiles`

#### Outputs
- A variety of plots and text files. These files will be located within the paths described in the `generalMutationAnalyses`  bucket from `config.json`

#### Execution

Linear execution:
```bash
cd generalMutationAnalyses
snakemake --keep-going --restart-times 2 --nolock 
```

Parallel execution in a SLURM cluster:
```bash
cd generalMutationAnalyses
``es 2 --keep-going --max-jobs-per-second 3 --max-status-checks-per-second 0.016 --cluster-config ../cluster.json --cluster-status jobState --jobs 500 --keep-going --cluster "../submit.py"
```

## Chromatin analyses

The following pipeline includes data analysis and plotting scripts for:
- Downloading signal chromatin maps for selected tissues and chromatin marks from the Roadmap epigenomics project.
- Processing of chromatin signal.
- Association analysis between chromatin and mutation load.

These analyses are included in Figure 2 of the publication.

##### Inputs
- Mutation maps and counts from [above](#mutation-calling)
- A variety of supporting files included in `auxiliaryFiles`

#### Outputs
- A variety of plots and text files. These files will be located within `$projectDir/chromatin/`, where `$projectDir` is the `projectDir` bucket in `config.json

#### Execution

Linear execution:
```bash
cd chromatin
snakemake --keep-going --restart-times 2 --nolock 
```

Parallel execution in a SLURM cluster:
```bash
cd chromatin
snakemake --restart-times 2 --keep-going --max-jobs-per-second 3 --max-status-checks-per-second 0.016 --cluster-config ../cluster.json --cluster-status jobState --jobs 500 --keep-going --cluster "../submit.py"
```

## Strand bias and selection analyses

The following pipeline includes data analysis and plotting scripts for:
- Processing and analysis of mutations in the transcribed and non-transcribed strands.
- Selection analyses using dN/dS.
- Selection analyses using VAF comparison.

These analyses are included in Figures 1,3 of the publication.

##### Inputs
- Mutation maps and counts from [above](#mutation-calling)
- A variety of supporting files included in `auxiliaryFiles`

#### Outputs
- A variety of plots and text files. These files will be located within the paths in the `selectionDir` bucket from `config.json`

#### Execution

Linear execution:
```bash
cd selectionAnalyses
snakemake --keep-going --restart-times 2 --nolock 
```

Parallel execution in a SLURM cluster:
```bash
``ctionAnalyses
snakemake --restart-times 1 --nolock --printshellcmds --keep-going --cluster-config ../cluster.json --cluster-status jobState --jobs 500 --cluster "../submit.py"
```

## Cancer analyses

The following pipeline includes data analysis and plotting scripts for:
- Enrichment of COSMIC mutations in mutation maps
- Mutations on cancer driver genes
- Selection analyses on cancer mutations and driver genes using dN/dS and VAF comparisons

These analyses are included in Figures 5 of the publication.

##### Inputs
- Mutation maps and counts from [above](#mutation-calling)
- A variety of supporting files included in `auxiliaryFiles`

#### Outputs
- A variety of plots and text files. These files will be located within the paths in the `cancer` bucket from `config.json`

#### Execution

Linear execution:
```bash
cd cancer
snakemake --keep-going --restart-times 2 --nolock 
```

Parallel execution in a SLURM cluster:
```bash
cd cancer
``atus-checks-per-second 0.016 --cluster-config ../cluster.json --cluster-status jobState --jobs 500 --keep-going --cluster "../submit.py"
```

