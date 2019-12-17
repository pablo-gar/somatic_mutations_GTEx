This repository contains all the software and scripts used for the following study.

**Garcia-Nieto PE, Morrison AJ, Fraser HB. The somatic mutation landscape of the human body. BioRxiv. 2019. [link](https://www.biorxiv.org/content/10.1101/668624v1)**

**Garcia-Nieto PE, Morrison AJ, Fraser HB. The somatic mutation landscape of the human body. Genome Biology. 2019. [link](https://www.biorxiv.org/content/10.1101/668624v1)**

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
    
***This document is not intended to be thouroguh description of methods or results. It is a guide that serves as a reproducibility reference for the code use in the aforementioned study.***

For a complete description of the methods and the results of theses analysis please refer to the publication

## Requirements
- Unix -- all software was run on bash
- snakemake -- for all pipelines
- STAR 2.5.2a (sequence aligner)
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

1. Setting up configuration files and downloading all required auxiliary files.
2. Creating a genome index for the STAR sequence aligner.
3. Mapping to the human genome.
4. Creation of readp pileup files for positions contating two different sequence calls.
5. Per-position addition of potential artifacts for the alternate allele.
6. Mutation calling and removal of likely false-positives.
7. Elimination of potential systematic artifacts (Panel of Normals).
8. Elimination of hyper-mutated samples.
9. High-level analyses

An in-depth description on how to execute each of these steps is described below.

## Seting up 

#### Downloading auxiliary files

Follow this [link](https://drive.google.com/a/stanford.edu/file/d/1v9ZIfkMmi7q8yh_lkn2BHFB4XDYxsERx/view?usp=sharing)

Uncompress it, and move it to a static location. These files are a variety of genome files, gene annotations, data from external databases (e.g. COSMIC, Roadmap epigenomics, TCGA), rna edit info, etc.

#### Configuring

Modify the file `config.json`, all buckets should be self explanatory. The file is setup as it was used for the original study, which was run on the Sherlock cluster at Stanford University.

It's worth noting the following:

- `genomeDir` is the path for the genome STAR index (see mapping for details in how to make it)
- Any bucket containing `\*/auxiliaryFiles/[...]` should be replace with the location of the downloaded and uncompress file above: `new\_path/auxiliaryFiles/[...]`
- `vcfFile` and `vcfFileCommon` are paths to the genotype vcf files from the GTEx project, please access those through dbGaP, feel free to contact me for any guidance on this.
- `projectDir` and `scratchDir` should be identical, this the root path where all results will be stored.
- Most other buckets can be left unmodified

#### Cluster and parallelization set up

Most pipelines are designed to run many jobs in parallel. Included there are some extra files and scripts to run snakemake on a SLURM cluster:


- The file 'cluster.json' contains rule-specific specifications and can be used on snakemake.
- submit.py is a custom submission script to be used for job submission on snakemake.
- jobState is a bash script to be used by snakemake to check the sate of jobs, **this script should be added to the bash `$PATH`**

## Creating genome index for STAR

This is the only step not included as a downdload or a step in any pipeline due to large memory and disk use.

#### Inputs
- Human genome Hg19 with common SNPs masked as "N". Included in the `auxiliaryFiles`, the path to the file should be the `genomeFasta` bucket in `config.json`
- Genome annotation in gft format. Included in the `auxiliaryFiles`, the path to the file should be the `auxiliaryFiles:GTF_gtex` bucket in `config.json`

#### Outputs
- Genome index for STAR

#### Execution
Replace variables ($) with the associated value in the `config.json` file

```bash
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ${genomeDir} --genomeFastaFile ${genomeFasta} --sjdbGTFfile ${auxiliaryFiles:GTF_gtex} --sjdbOverhang 99
```

## Mapping raw RNA-seq reads, and creating pileup files

The following pipeline includes:
- Mapping to the human genome, sorting and indexing bam files, removing duplicate reads.
- Creating coverage maps.
- Creating read pipelup files for positions covered by two difference sequence calls
- Adding per-site bias information to pileup files: read position bias, mapping quality bias, sequence quality bias, strand bias, variant distance bias. 

A full description of these biases can be found in the Methods section of the publication.

#### Inputs
- Raw rna-seq paired fastq files. Included in the `auxiliaryFiles` download, the path to the file should be the `fastqDir` bucket in `config.json`. File names should have this format `{sample}_1.fastq.gz` `{sample}_2.fastq.gz`
- SRA table with sample information. Included in the `auxiliaryFiles` download, the path to the file should be the `sraTable` bucket in `config.json`
- Genome size table. Included in the `auxiliaryFiles` download, the path to the file should be the `genomeSizeFile` bucket in `config.json`
- A variety supporting files included in `auxiliaryFiles`

#### Outputs
- Sorted, without duplicate reads bam files. Files will be locate in `$mappingDir/{tissue}/{sample}_RmdupSortedAligned.out.bam`, where `$mappingDir` is the `mappingDir` bucket in `config.json`
- Coverage maps `$root/depth_bam/{sample}.bed.gzip`. Where `$root` is the `projectDir` bucket in `config.json`
- Pileup files for positions with two sequence calls `$pileupDir/{tissue}/{sample}.txt`. Where `$pileupDir` is relative and is the `pileupDir` bucket in `config.json`

#### Execution

Linear execution (this is not feasible as it would take a very large time to finish):
```bash
cd mappingMutationCalling_pipeline
snakemake --snakefile SnakefilePileup.smk --printshellcmds --keep-going --restart-times 2 
```

Parallel execution:
```bash
cd mappingMutationCalling_pipeline
snakemake --snakefile SnakefilePileup.smk --printshellcmds --keep-going --max-jobs-per-second 3 --max-status-checks-per-second 0.016 --nolock --restart-times 2 --cluster-config ../cluster.json --cluster-status jobState --jobs 500 --cluster "../submit.py"
```

## Mutation calling

The following pipeline includes:
- Calling mutations based on coverage and number of read supporting alt allele, while incorporating sequencing error probabilties.
- Label and remove false positives:
    - Black listed regions.
    - RNA edits.
    - Splicing junction errors.
    - Sequencing errors.
    - All biases described in previous section 
- Remove mutations from systematic artifacts based on a pseudo Panel of Normals (PoN)
- Remove hypermutated samples

A full description of the false-positive removal can be found in the Methods section of the publication.

#### Inputs
- Pileup files for positions with two sequence calls `$pileupDir/{tissue}/{sample}.txt`. *From previous section*
- A variety supporting files included in `auxiliaryFiles`

#### Outputs
- Mutation maps per sample `${mutationCountDir:map}/{tissue}/n6_0.0_0.7/{sample}.txt`.
- Counts of mutations following a 6-type profile (C\>T, C\>G, C\>T, C\>A, T\>G, T\>A, T\>C) `$mutationCountDir:$context/{tissue}/n6_0.0_0.7/{sample}.txt`
- Counts of mutations with pentanucleotide context per sample `$mutationCountDir:$context/{tissue}/n6_0.0_0.7/{sample}.txt`

#### Execution

Linear execution (this is not feasible as it would take a very large time to finish):
```bash
cd mutationCalling
snakemake --keep-going --restart-times 2 --nolock 
```

Parallel execution:
```bash
cd mutationCalling
snakemake --keep-going --max-jobs-per-second 15 --restart-times 2 --max-status-checks-per-second 0.016 --nolock --cluster-config ../cluster.json --cluster-status jobState --jobs 500 --cluster "../submit.py"
```
