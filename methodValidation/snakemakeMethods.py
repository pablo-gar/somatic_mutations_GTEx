import os, subprocess
from time import gmtime, strftime

def getPartionedGenome(genomeSizeFile, regionSize):
        
    ''' 
    Returns a list of list, with inner list containing regions of the genome with the format
    [chr, start, end]
    Upper list ranges from 1 to number of jobs per bam
    Second lists are regions homogeneously distributed
    Third list contains the actual region info in the format indicated above
    '''
    
    regions = []
    regionSize = int(regionSize)

    genomeSize = open(genomeSizeFile, 'r')
    for chrom in genomeSize:
        chrom = chrom.rstrip().split("\t")
        nwindows = round((float(chrom[1]) / regionSize) + 0.5)
        start = 1 
        for i in range(nwindows):
            end = start + regionSize
            if end > int(chrom[1]): end = int(chrom[1]) - 1 
            regions.append([chrom[0], str(start), str(end)])
            start = start + regionSize + 1 
            #if(i > 3):
            #    break
        #break
    
    genomeSize.close()
    
    return regions

def repeatList(x, n):
    
    ''' Repeats n times all elements of list x'''
    
    return [repeated for repeated in x for i in range(n)]

def expandPaths(obj):
    
    '''
    Expands all paths in the leaves of a dictionary of dictionaries (i.e. config object in snakemake
    nodes of tree have to be dictionaries and leaves have to be strings
    '''
    
    if isinstance(obj, str):
        obj = os.path.expanduser(obj)
    else:
        for i in obj:
            obj[i] = expandPaths(obj[i])
    
    return(obj)

def getSampleInVcf(samples, vcf):
    
    '''
    Returns a boolean list with true in position corresponding to samples
    present in the vcf
    '''
    cmd = "bcftools query -l " + vcf
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    samplesVcf = process.communicate()[0].decode("utf-8").split("\n")
    result = []
    for s in samples:
        result.append(s in samplesVcf)
    
    return result
    
def sraToGtex(sra, sraTable):
    cmd = "grep -E " + sra + " " + sraTable + " | cut -f 32"
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    return(process.communicate()[0].decode("utf-8").rstrip())

def writeTimeToFile(timeFile):
    TIME =  strftime("%Y-%m-%d_%H:%M:%S.", gmtime())
    lastTime = open(timeFile, "w")
    lastTime.write(TIME)
    lastTime.close()
    return(TIME)
    
def readTimeFromFile(timeFile):
    lastTime = open(timeFile, "r")
    TIME = lastTime.readline()
    lastTime.close()
    return TIME
