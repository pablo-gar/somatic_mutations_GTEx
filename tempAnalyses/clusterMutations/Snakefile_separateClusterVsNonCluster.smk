# Run line

import sys, os


##--------------------------------------------------##
## GLOBAL

configfile: "../../config.json"

ROOT_DIR = "/scratch/users/paedugar/somaticMutationsProject"
WORKING_DIR = os.path.join(ROOT_DIR, "clusterMutations")

MAF = "n6_0.0_0.5"

# Sample ids
TISSUES, SAMPLES, = glob_wildcards(os.path.join(ROOT_DIR, config["mutationCountDir"]["root"], "filtered", '{tissue,\w+}', MAF, '{sample,\w+}.txt'))
#TISSUES = TISSUES[0:2]
#SAMPLES = SAMPLES[0:2]


# COMMENT THIS if you get tissues and MAF from folder names 
MAF = [MAF] * len(SAMPLES)
# FINISH COMMENTING here

clusterTypes = ["cluster", "nonCluster"]


##--------------------------------------------------##


#-------------------------------------------
# RULEs to separate cluster vs non-clustered mutations
# Once mutations have been called and a filter column has been appended
rule all:
    input:
        # Cluster
        expand(os.path.join(WORKING_DIR, "cluster", config["mutationCountDir"]["count"], '{tissue}', '{maf}', '{sample}.txt'), zip, tissue = TISSUES, maf = MAF, sample =  SAMPLES),
        expand(os.path.join(WORKING_DIR, "cluster", config["mutationCountDir"]["context"], '{tissue}', '{maf}', '{sample}.txt'), zip, tissue = TISSUES, maf = MAF, sample =  SAMPLES), 
        expand(os.path.join(WORKING_DIR, "cluster", config["mutationCountDir"]["map"], '{tissue}', '{maf}', '{sample}.txt'), zip, tissue = TISSUES, maf = MAF, sample =  SAMPLES),
        # Non cluster
        expand(os.path.join(WORKING_DIR, "nonCluster", config["mutationCountDir"]["count"], '{tissue}', '{maf}', '{sample}.txt'), zip, tissue = TISSUES, maf = MAF, sample =  SAMPLES),
        expand(os.path.join(WORKING_DIR, "nonCluster", config["mutationCountDir"]["context"], '{tissue}', '{maf}', '{sample}.txt'), zip, tissue = TISSUES, maf = MAF, sample =  SAMPLES), 
        expand(os.path.join(WORKING_DIR, "nonCluster", config["mutationCountDir"]["map"], '{tissue}', '{maf}', '{sample}.txt'), zip, tissue = TISSUES, maf = MAF, sample =  SAMPLES)

rule create_clusterFiles:
    input:
        os.path.join(ROOT_DIR, config["mutationCountDir"]["root"], "filtered", '{tissue}', '{maf}', '{sample}.txt')
    params:
        contextLength = 5
    output:
        count = os.path.join(WORKING_DIR, "cluster", config["mutationCountDir"]["count"], '{tissue}', '{maf}', '{sample}.txt'),
        context = os.path.join(WORKING_DIR, "cluster", config["mutationCountDir"]["context"], '{tissue}', '{maf}', '{sample}.txt'),
        mapOut = os.path.join(WORKING_DIR, "cluster", config["mutationCountDir"]["map"], '{tissue}', '{maf}', '{sample}.txt')
    shell:
        '''
        # Filter
        awk -v FS="\\t" -v OFS="\\t" '{{if ($8 == "clustered_mutation") print $0}}' {input} | cut -f-7 > {output.mapOut}
        
        module load R/3.4.0
        Rscript ../../mutationCalling/R/printCountsContext_fromMap.R {params} {output.mapOut} {output.count} {output.context}
        '''
        
rule create_NonClusterFiles:
    input:
        os.path.join(ROOT_DIR, config["mutationCountDir"]["root"], "filtered", '{tissue}', '{maf}', '{sample}.txt')
    params:
        contextLength = 5
    output:
        count = os.path.join(WORKING_DIR, "nonCluster", config["mutationCountDir"]["count"], '{tissue}', '{maf}', '{sample}.txt'),
        context = os.path.join(WORKING_DIR, "nonCluster", config["mutationCountDir"]["context"], '{tissue}', '{maf}', '{sample}.txt'),
        mapOut = os.path.join(WORKING_DIR, "nonCluster", config["mutationCountDir"]["map"], '{tissue}', '{maf}', '{sample}.txt')
    shell:
        '''
        # Filter
        awk -v FS="\\t" -v OFS="\\t" '{{if ($8 == "PASS") print $0}}' {input} | cut -f-7 > {output.mapOut}
        
        module load R/3.4.0
        Rscript ../../mutationCalling/R/printCountsContext_fromMap.R {params} {output.mapOut} {output.count} {output.context}
        '''
