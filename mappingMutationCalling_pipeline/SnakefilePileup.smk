
import os, pandas, re, sys
import snakemakeMethods as SM
from numpy import repeat

##--------------------------------------------------##
## GLOBAL

configfile: "../config.json"
localrules: check_mutationPileup, organize_tissues, all
ruleorder: merge_pileups > create_mutationPileup

# Expand ~ to /home/user, or whatever the system default is
config = SM.expandPaths(config)
 
# Looks for all the fastq files in the fastqDir and subfolders
CURRENT_SAMPLES = [re.sub(r"_\d+.fastq.gz", "", x) for root, dirs, files in os.walk(config["fastqDir"]) for x in files ]
CURRENT_SAMPLES = set(CURRENT_SAMPLES)

# Gets what tissue the samples are form
SAMPLES = pandas.read_table(config["sraTable"], index_col = 15)[['body_site_s', 'submitted_subject_id_s']]
SAMPLES = SAMPLES.assign(Run_s = SAMPLES.index.values)
SAMPLES['body_site_s'] = [re.sub(r"\s+", "_", re.sub(r"-", " ", re.sub(r"[\(\)]", "", x))) for x in SAMPLES['body_site_s']]
SAMPLES = SAMPLES[ [x in CURRENT_SAMPLES for x in SAMPLES.Run_s] ]

# Genome partioned
GENOME_PARTITIONS =  [".".join(x) for x in SM.getPartionedGenome(config["genomeSizeFile"], config["pileupPartitionSize"])]



# Suffix for corresponding paired-read files
PAIR1 = '{sample}_1.fastq.gz'
PAIR2 = '{sample}_2.fastq.gz'

##--------------------------------------------------##


##--------------------------------------------------##
## PIPELINE

rule all:
    input: 
        # Create joint pileup files
        expand(os.path.join(config["projectDir"], config["pileupDir"], '{tissue}', '{sample}.txt'), zip, tissue = SAMPLES.body_site_s, sample = SAMPLES.Run_s),
    

rule star_map:
    input:
        p1 = os.path.join(config["fastqDir"], "{tissue}", PAIR1),
        p2 = os.path.join(config["fastqDir"], "{tissue}", PAIR2),
        genomeDir = config["genomeDir"]
    params:
        miss = config["missMatchesMapping"]
    threads: 8
    output: 
        temp(os.path.join(config["mappingDir"], "{tissue}", '{sample}_Aligned.out.bam'))
    shell:
        """
            rm -rf {config[scratchDir]}/{config[tempDir]}/{wildcards.sample}
            STAR \
            --runThreadN {threads} \
            --genomeDir {input.genomeDir}\
            --readFilesCommand zcat \
            --outFilterMultimapNmax 1 \
            --outSAMtype BAM Unsorted \
            --outSAMattributes MD NH \
            --clip5pNbases 6 \
            --outFilterMismatchNmax {params.miss} \
            --readFilesIn {input.p1} {input.p2} \
            --outTmpDir {config[scratchDir]}/{config[tempDir]}/{wildcards.sample} \
            --outFileNamePrefix {config[mappingDir]}/{wildcards.tissue}/{wildcards.sample}_
        """

rule sort_bam:
    input:
        os.path.join(config["mappingDir"], '{tissue}', '{sample}_Aligned.out.bam')
    output:
        temp(os.path.join(config["mappingDir"], '{tissue}', '{sample}_SortedAligned.out.bam'))
    threads: 1
    shell:
        """
            rm -rf {config[mappingDir]}/{wildcards.sample}*
        
            {config[samtoolsBin]} \
            sort \
            -l 0 \
            -@ {threads} \
            -T {config[mappingDir]}/{wildcards.sample} \
            -O BAM \
            -l 9 \
            -o {output} \
            {input}
        """

rule removeDuplicates_bam:
    input:
        os.path.join(config["mappingDir"], '{tissue}', '{sample}_SortedAligned.out.bam')
    output:
        os.path.join(config["mappingDir"], '{tissue}', '{sample}_RmdupSortedAligned.out.bam')
    threads: 1
    shell:
        'python/removeDuplicatesBam.py -i {input} -o {output}'

rule index_bam:
    input:
        os.path.join(config["mappingDir"], '{tissue}', '{sample}_RmdupSortedAligned.out.bam')
    output:
        os.path.join(config["mappingDir"], '{tissue}', '{sample}_RmdupSortedAligned.out.bam.bai')
    threads: 1
    shell:
        '{config[samtoolsBin]} index {input}'

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

rule create_mutationPileup:
    input: 
        bam = os.path.join(config["mappingDir"], '{tissue}', '{sample}_RmdupSortedAligned.out.bam'),
        bamIndex = os.path.join(config["mappingDir"], '{tissue}', '{sample}_RmdupSortedAligned.out.bam.bai'),
        coverageMap = os.path.join(config["scratchDir"], 'depth_bam', '{tissue}', '{sample}.bed.gzip'), # This is here just to make sure that the coverage maps are done
        genomeFasta = config["genomeFasta"]
         
    params:
        # For calling mutations
        chromosome = lambda wildcards: wildcards.region.split(".")[0],
        start = lambda wildcards: wildcards.region.split(".")[1],
        end = lambda wildcards: wildcards.region.split(".")[2], # This is here just to make sure that the coverage maps are done
        baseQualityCutoff = config["seqQualityCutoff"],
        commonSNPs = "N",
        checkGenotype = "N",
        out_with_germline = os.path.join(config["projectDir"], config["pileupDirWithGermline"], '{tissue}', '{sample}_{region}.txt'),
        # For eliminating germline
        sampleId = lambda wildcards: SAMPLES["submitted_subject_id_s"][wildcards.sample],
        vcfFile = config["vcfFile"],
        tempDir = os.path.join(config["scratchDir"], config["tempDir"])
         
    output:
        os.path.join(config["projectDir"], config["pileupDir"], '{tissue}', '{sample}_{region}.txt')
         
    shell:
        '''
        
        # Create pileup
        python3 python/SNVinBam_pileup.py {input.genomeFasta} {input.bam} {params.chromosome} {params.start} {params.end} {params.baseQualityCutoff} {params.commonSNPs} {params.checkGenotype} | awk -v FS="\\t" -v OFS="\\t" '{{$2=$2 "\\t" $2; print $0}}' | sort -k1,1 -k2,2n > {params.tempDir}/{wildcards.sample}_{wildcards.region}.no_filters

        
        # Append false-positive probabilities if file is not empty
        mkdir -p $(dirname {params.out_with_germline})
        
        if [ -s  {params.tempDir}/{wildcards.sample}_{wildcards.region}.no_filters ]
        then
        
            bcftools mpileup -B --skip-indels --min-BQ {params.baseQualityCutoff} --max-depth 50000 --ignore-RG -T {params.tempDir}/{wildcards.sample}_{wildcards.region}.no_filters --fasta-ref {config[genomeFasta]} {input.bam} > {params.tempDir}/{wildcards.sample}_{wildcards.region}.vcf
            
            bcftools query -f '%CHROM\\t%POS\\t%POS\\t[%INFO/VDB]\\t[%INFO/RPB]\\t[%INFO/MQB]\\t[%INFO/BQB]\\t[%INFO/MQSB]\\n' {params.tempDir}/{wildcards.sample}_{wildcards.region}.vcf | sort -k1,1 -k2,2n > {params.tempDir}/{wildcards.sample}_{wildcards.region}.stats
            
            rm {params.tempDir}/{wildcards.sample}_{wildcards.region}.vcf
            
            bedtools intersect -sorted -f 1 -wa -wb -a {params.tempDir}/{wildcards.sample}_{wildcards.region}.no_filters -b {params.tempDir}/{wildcards.sample}_{wildcards.region}.stats | cut -f3,12,13,14 --complement > {params.out_with_germline}
            
            rm {params.tempDir}/{wildcards.sample}_{wildcards.region}.stats
            rm {params.tempDir}/{wildcards.sample}_{wildcards.region}.no_filters
            
        else
        
            mv  {params.tempDir}/{wildcards.sample}_{wildcards.region}.no_filters {params.out_with_germline}
            
        fi
        
        # Remove germline
        python3 python/SNV_removeGermiline_correctContext.py {params.out_with_germline} {params.sampleId} {params.vcfFile} {params.tempDir} > {output}
        rm {params.out_with_germline}
        '''

rule check_mutationPileup:
    input:
        expand(os.path.join(config["projectDir"], config["pileupDir"], '{tissue}', '{sample}_{region}.txt'), zip, \
                                                                                       tissue = SM.repeatList(SAMPLES.body_site_s, len(GENOME_PARTITIONS)), \
                                                                                       sample = SM.repeatList(SAMPLES.Run_s, len(GENOME_PARTITIONS)), \
                                                                                       region = GENOME_PARTITIONS * len(SAMPLES) \
              )
    output:
        temp(os.path.join(config["scratchDir"], config["tempDir"], 'done_mutationPileup.txt'))
    shell:
        'echo DONE > {output}'
        
rule merge_pileups:
    input:
        ancient(os.path.join(config["scratchDir"], config["tempDir"], 'done_mutationPileup.txt'))
    output:
        os.path.join(config["projectDir"], config["pileupDir"], '{tissue}', '{sample}.txt')
    shell: 
        '''
        allRegions={config[projectDir]}/{config[pileupDir]}{wildcards[tissue]}/{wildcards[sample]}_*.txt
        if ls $allRegions 1> /dev/null 2>&1
        then 
            cat $allRegions | sort -k1,1 -k2,2n > {output} && rm $allRegions
        fi
        '''
       
