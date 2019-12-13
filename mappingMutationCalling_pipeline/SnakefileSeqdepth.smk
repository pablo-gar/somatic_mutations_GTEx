
import os, sys
import snakemakeMethods as SM

##--------------------------------------------------##
## GLOBAL

configfile: "../config.json"

# Expand ~ to /home/user, or whatever the system default is
config = SM.expandPaths(config)
 
# 
TISSUES, SAMPLES, = glob_wildcards(os.path.join(config["mappingDir"], '{tissue,\w+}', '{sample,\w+}_RmdupSortedAligned.out.bam'))
##--------------------------------------------------##


##--------------------------------------------------##
## PIPELINE

rule all:
    input: 
        expand(os.path.join(config["scratchDir"], 'depth_bam', '{tissue}', '{sample}.bed.gzip'), zip, tissue = TISSUES, sample = SAMPLES)

rule create_coverageMaps:
    input:
        os.path.join(config["mappingDir"], '{tissue}', '{sample}_RmdupSortedAligned.out.bam')
    params:
        seqQual = int(config["seqQualityCutoff"]) - 1,
        coverage = config["coverageCutoff"],
        tempDir = os.path.join(config["scratchDir"], config["tempDir"])
    output:
        os.path.join(config["scratchDir"], 'depth_bam', '{tissue}', '{sample}.bed.gzip')
    shell:
        '''
        samtools depth -q {params.seqQual} -d 0 {input} | awk -v FS="\\t" -v OFS="\\t" '{{if ($3 >= {params.coverage}) print $1 "\\t" $2 "\\t" $2 + 1 "\\tNA\\t" $3 "\\t*" }}' | sort -k 1,1 -k2,2n | gzip --best -c > {output}
        '''
