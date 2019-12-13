import os
import snakemakeMethods as SM
##--------------------------------------------------##
## GLOBAL
#
localrules: all
configfile: "../config.json"
config = SM.expandPaths(config)


DIR_NAME = "method_validation"
MUT_DIR = os.path.join(DIR_NAME, "mutations")
EXO_FIRST = os.path.join(MUT_DIR, "exome_first")
RNA_FIRST = os.path.join(MUT_DIR, "RNA_first")
RNA_FIRST_context = os.path.join(MUT_DIR, "RNA_first_context")

MUT_DIR_LAX = os.path.join(DIR_NAME, "mutations_lax")
EXO_FIRST_LAX = os.path.join(MUT_DIR_LAX, "exome_first")
RNA_FIRST_LAX = os.path.join(MUT_DIR_LAX, "RNA_first")
RNA_FIRST_LAX_context = os.path.join(MUT_DIR_LAX, "RNA_first_context")

FDR = os.path.join(DIR_NAME, "fdr")

DEPTH_DIR = "depth_bam"

FILES, = glob_wildcards(os.path.join(config["projectDir"], RNA_FIRST_LAX, '{FILES,.+}.txt'))

# Reads conversion table for sample ids



##--------------------------------------------------##


##--------------------------------------------------##
## PIPELINE

rule all:
    input: 
        ##Lax mutations with contex
        expand(os.path.join(config["projectDir"], RNA_FIRST_LAX_context, '{file}.txt'), file = FILES)
        

rule get_context:
    input:
        matched = os.path.join(config["projectDir"], RNA_FIRST_LAX, '{file}.txt'),
        genome_bedtools = os.path.join(config["auxiliaryFiles"]["genomeFolder"], "GRCh37.p13.genome.encode.genomeFile.bedtools"),
        genome_fasta =  os.path.join(config["auxiliaryFiles"]["genomeFolder"], "GRCh37.p13.genome.encode.fa"),
        exon_bed =  os.path.join(config["auxiliaryFiles"]["dir"], "gencode.v19.genes.v7.patched_contigs.onlyExons.bed")
    params:
        sample = "{file}",
        temp_dir =  os.path.join(config["projectDir"] ,config["tempDir"])
    output:
        os.path.join(config["projectDir"], RNA_FIRST_LAX_context, '{file}.txt')
    shell:
        """
        grep $'\\t+\\t' {input.exon_bed} | \
            bedtools intersect -f 1 -wa -a <(awk -v OFS="\\t" '$2 = $2-1 "\\t" $2' {input.matched}) -b stdin | \
            bedtools slop -b 1 -g {input.genome_bedtools}  -i stdin | \
            bedtools getfasta -fi {input.genome_fasta} -bedOut -bed stdin | \
            awk -v OFS="\\t" '{{print $0 "\\t+"}}' > {params.temp_dir}/{params.sample}_1.txt

        grep $'\\t-\\t' {input.exon_bed} | \
            bedtools intersect -f 1 -wa -a <(awk -v OFS="\\t" '$2 = $2-1 "\\t" $2' {input.matched}) -b stdin | \
            bedtools slop -b 1 -g {input.genome_bedtools}  -i stdin | \
            bedtools getfasta -fi {input.genome_fasta} -bedOut -bed stdin | \
            awk -v OFS="\\t" '{{print $0 "\\t-"}}' > {params.temp_dir}/{params.sample}_2.txt

        cat {params.temp_dir}/{params.sample}_1.txt {params.temp_dir}/{params.sample}_2.txt | sort -k1,1 -k2,2n | awk -v OFS="\\t" '$2=$2+2' | cut -f1-2,4- > {output}
        """ 
