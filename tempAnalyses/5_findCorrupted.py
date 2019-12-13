#!/usr/bin/env python3
# Scans SRA fastq files from the master table and finds out if the number 
# of lines is consistent
# If it is not, it deletes the files
#
# For all files creates a table and annotates if it is consistent or not

#python3 5_findCorrupted.py /scratch/PI/hbfraser/gtex/anno/GtexSraRunTable.txt /scratch/users/paedugar/somaticMutationsProject/fastq /scratch/users/paedugar/somaticMutationsProject/auxiliaryFiles/integrityFastq.txt

import pandas, sys, os, re, fyrd

ASSAYS = {"Allele-Specific Expression":"RNASEQ", "RNA Seq (NGS)":"RNASEQ", "Whole Genome (NGS)":"WGS", "Whole Exome (NGS)":"EXO", "SNP/CNV Genotypes (NGS)":"EXO", "CNV Genotypes":"EXO"}


def main():
    
    sampleTableFile = sys.argv[1]
    dataFolder = sys.argv[2]
    outTableFile = sys.argv[3]
    outDir = os.path.dirname(outTableFile)
    
    
    # Creating output table if does not exits
    if not os.path.exists(outTableFile):
        outTable = open(outTableFile, "w")
        outTable.write("\t".join(["Run_s", "Corrupted", "Reason"]) + "\n")
        outTable.close()
        
    outTable = pandas.read_table(outTableFile, sep = "\t", header = 0)
    outTable = outTable.set_index('Run_s', drop = False)
    
    # Reading table of samples
    sampleTable = pandas.read_table(sampleTableFile, index_col = 15, sep = "\t")[['body_site_s', 'molecular_data_type_s' ]]
    
    # Looping of all samples and looking for fastq files
    jobs = []
    for index, row in sampleTable.iterrows():
        
        tissue = re.sub(r"\s+", "_", re.sub(r"-", " ", re.sub(r"[\(\)]", "",  row['body_site_s'])))
        experiment =  ASSAYS[row['molecular_data_type_s']]
        sraId = index
        samplePath_read1 = os.path.join(dataFolder, tissue, "_".join([sraId, "1.fastq.gz"]))
        samplePath_read2 = os.path.join(dataFolder, tissue, "_".join([sraId, "2.fastq.gz"]))
        
            
        # If fastq file exists and it has not been checked, then check for integreity
        if os.path.exists(samplePath_read1) and os.path.exists(samplePath_read2):
            
            # If it has not been checked before
            if not sraId in outTable.index.values:
                jobs.append(checkForCorruption(sraId, samplePath_read1, samplePath_read2, os.path.join(outDir, "integrity_intermediate_" + sraId)))
            
                
            # If it was checked before but corrupted check again
            if sraId in outTable.index.values and outTable["Corrupted"][sraId]:
                jobs.append(checkForCorruption(sraId, samplePath_read1, samplePath_read2, os.path.join(outDir, "integrity_intermediate_" + sraId)))
                
    # Waiting for jobs to finish
    for job in jobs:
        job.wait()
        job.clean()
    
    # Putting intermidate tables together
    os.system("cat "+ outTableFile + " " + os.path.join(outDir, "integrity_intermediate_*") + " > " + os.path.join(outDir, "temp"))
    os.rename(os.path.join(outDir, "temp"), outTableFile)
    os.system("rm " + os.path.join(outDir, "integrity_intermediate_*"))
    

def checkForCorruption(sraId, samplePath_read1, samplePath_read2, outTableFile):
    command = ["module load fraserconda", 
               " ".join(["/home/users/paedugar/scripts/FraserLab/somaticMutationsProject/tempAnalyses/5_findCorrupted_Job.py", outTableFile, sraId, samplePath_read1, samplePath_read2])
               ]
    command = "\n".join(command)
    job = fyrd.Job(command, nodes = 1, cores = 8, partition = "hbfraser,hns,normal,owners")
    job.submit()
    return job
    
    
if __name__=="__main__":
    main()
