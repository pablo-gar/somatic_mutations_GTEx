
import os, pandas, re
import snakemakeMethods as SM
from numpy import repeat

##--------------------------------------------------##
## GLOBAL

configfile: "../config.json"
localrules: organize_tissues, all

# Expand ~ to /home/user, or whatever the system default is
config = SM.expandPaths(config)
 
# Looks for all the fastq files in the fastqDir and subfolders
CURRENT_SAMPLES = set([re.sub(r"_\d+.fastq.gz", "", x) for root, dirs, files in os.walk(config["fastqDir"]) for x in files ])

# Gets what tissue the samples are form
SAMPLES = pandas.read_table(config["sraTable"], index_col = 15)[['body_site_s', 'submitted_subject_id_s']]
SAMPLES = SAMPLES.assign(Run_s = SAMPLES.index.values)
SAMPLES['body_site_s'] = [re.sub(r"\s+", "_", re.sub(r"-", " ", re.sub(r"[\(\)]", "", x))) for x in SAMPLES['body_site_s']]
SAMPLES = SAMPLES[ [x in CURRENT_SAMPLES for x in SAMPLES.Run_s] ]

# Genome partioned
GENOME_PARTITIONS =  [".".join(x) for x in SM.getPartionedGenome(config["genomeSizeFile"], config["pileupPartitionSize"])]
region = GENOME_PARTITIONS


# Suffix for corresponding paired-read files
PAIR1 = '{sample}_1.fastq.gz'
PAIR2 = '{sample}_2.fastq.gz'

##--------------------------------------------------##


##--------------------------------------------------##
## PIPELINE

rule all:
    input: 
       expand(os.path.join(config["fastqIntegrityDir"], "{tissue}", "{sample}.txt"), zip, tissue=SAMPLES.body_site_s, sample = SAMPLES.Run_s)

        
rule organize_tissues:
    input:
        p1 = os.path.join(config["fastqDir"],  PAIR1),
        p2 = os.path.join(config["fastqDir"],  PAIR2)
    output:
        p1 = os.path.join(config["fastqDir"], "{tissue}", PAIR1),
        p2 = os.path.join(config["fastqDir"], "{tissue}", PAIR2)
    shell:
        "mv {input.p1} {output.p1} && mv {input.p2} {output.p2}"
        
rule check_fastqIntegrity:
    input: 
        p1 = os.path.join(config["fastqDir"],  "{tissue}", PAIR1),
        p2 = os.path.join(config["fastqDir"],  "{tissue}", PAIR2)
    params:
        lambda wildcards: os.path.join(config["fastqPurgatoryDir"], wildcards.tissue)
    output:
        os.path.join(config["fastqIntegrityDir"], "{tissue}", "{sample}.txt")
    shell:
        '''
        size1=$(pigz -d -c -p 1 {input.p1} |  wc -l)
        size2=$(pigz -d -c -p 1 {input.p2} |  wc -l)

        if [ $size1 == $size2 ]
        then
            echo -e "correct\t$size1\t$size2" > {output}
        else
            mkdir -p {params}
            mv {input.p1} {params}/
            mv {input.p2} {params}/
            echo -e "corrupted\t$size1\t$size2" > {output}
        fi
        '''
