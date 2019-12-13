#Takes a pileup file and removes the positions that fall in a germiline
#variants according to a vcfFile. It also corrects whenever there are Ns
#in the context
#
#Author: Pablo E. Garcia-Nieto
#Date: 01/03/2018
#
#Requires an executable of ov bcftools
#
# @args[1] pileup file
# @args[2] smaple Id
# @args[3] vcfFile
# @args[4] folder to stores temp files
#
# Prints result to STDOUT

import sys, os, re, random, string, subprocess

CONTEXT_I = 9

def main():
    
    pileupFile, sampleId, vcfFile, tempFolder = readArgs()
    
    contextLength = getContextLength(pileupFile)
    
    # Creates a file that contains all the positions to look for SNPs
    positionsFile = createPositionsFile(pileupFile, contextLength, tempFolder)
    
    # Creates a matrix with containing all SNPs
    genotypeDict = createGenotypeFile(positionsFile, sampleId, vcfFile)
    os.remove(positionsFile)
    
    # Goes over pileup and ignores positions with germiline as well as corrects context when possible
    printCorrected(pileupFile, contextLength, genotypeDict)
    
    
def printCorrected(pileupFile, contextLength, genotypeDict):
    
    pileup = open(pileupFile, "r")
    
    for line in pileup:
        
        line = line.rstrip().split("\t")
        context = line[CONTEXT_I]
        chrom = line[0]
        pos = line [1]
        
        if assertGenotype(genotypeDict, chrom, pos):
            #Correct context if possible
            line[CONTEXT_I] = correctContext(line[CONTEXT_I], genotypeDict, chrom, int(pos) - round(contextLength / 2 - 0.5))
            line[2] = line[CONTEXT_I][round(contextLength / 2 - 0.5)]
            print(*line, sep = "\t")
        
def correctContext(refBase, genotypeDict, chrom, refPos):

    '''
    Given a genomic sequence and a vcfFile, it will get assigned the correct genotype
    where possible in the Ns
    
    refPos is 0-based
    vcfFile is 1-based
    
    @return corrected sequence
    '''
    
    chrom = re.sub("chr", "", chrom)
    
    # Get all genoypes in this region
    genotypes = set()
    for i in range(refPos, refPos + len(refBase)):
        region = ".".join([chrom, str(i)])
        if region in genotypeDict:
            genotypes = genotypes.union(genotypeDict[region])

    while(len(genotypes) > 0):
        genotype = genotypes.pop()
        old = genotype
        # Parses out genotype info out of bcftools output
        genotype = genotype.split("|")
        pos = int(genotype[1])
        refAlt = genotype[2:4] #Alleles
        genotype = genotype[4].split("/")
        
        if len(refAlt[0]) > 1 or len(refAlt[1]) > 1: continue

        if genotype[0] != genotype[1]: continue
        if genotype[0] not in ["0","1"]: continue

        genotype = [int(i) for i in genotype]
        allele =  refAlt[genotype[0]] # Only one allele
        
        c = pos - refPos
        refBase = "".join([refBase[0:c], allele, refBase[c + 1:]])

    return refBase
 
def assertGenotype(genotypeDict, chrom, pos):
    ''' 
    Asserts whether this position passes a genotype quality controls
    Essentially asserts whether there is a heterozygous germline variant
    
    Returns True if:
        - This position does not have an annotation in the vcf file
        - This position does not have a het genotype for the sample
    TODO
        - This position contains at least one read for the two alleles for the genotype
    '''
    
    assertion = True
    chrom = re.sub("chr", "", chrom)
    region = ".".join([chrom, pos])
    
    if region in genotypeDict:
        
        genotypes = [x.split("|")[4] for x in genotypeDict[region]] # Get genotypes 
        for i in genotypes:
            if i == "0/1" or i == "1/0" or i == "1/." or i == "0/." or i == "./1" or i == "./0":
                assertion = False
    
    return assertion
    
    
def readArgs():
    
    return sys.argv[1:5]

def getContextLength(pileupFile):
    
    with open(pileupFile, "r") as pileup:
        line = pileup.readline()
        if line == "": sys.exit()
        line = line.rstrip().split("\t")
        context = line[CONTEXT_I]
        
        return len(context)

def createPositionsFile(pileupFile, contextLength, tempFolder):
    
    pileup = open(pileupFile, "r")
    tempPosFile = os.path.join(tempFolder, os.path.basename(pileupFile) + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(20)))
            
    tempPos = open(tempPosFile, "w")
    
    for line in pileup:
        line = line.rstrip().split("\t")
        chrom = re.sub(r"chr", "", line[0])
        pos = line[1]
        posStart = str(int(pos) - round(contextLength / 2 - 0.5))
        posEnd = str(int(pos) + round(contextLength / 2 - 0.5))
        
        tempPos.write("\t".join([chrom, posStart, posEnd, "\n"]))
        
    pileup.close()
    tempPos.close()
    
    return(tempPosFile)


def createGenotypeFile(positionsFile, sampleId, vcfFile):
    
    '''
    Saves output in dictionary with keys being chr.pos
    contents of dict is a set with each element containing a genotype line for that pos
    '''
    
    
    bcfToolsCmd = ["bcftools", "query", "-f", "%CHROM|%POS|%REF|%ALT|[%GT]\n", "-s", sampleId, "-R", positionsFile, vcfFile]
    bcfTools = subprocess.Popen(bcfToolsCmd, stdout = subprocess.PIPE)
    
    genotype = dict()
    while True:
        output = bcfTools.stdout.readline().decode("ascii")
        if output == '' and bcfTools.poll() is not None:
            break
        if output:
            output = output.rstrip()
            chrom, originalPos, refAllele = output.split("|")[0:3]
            
            
            for i in range(len(refAllele)):
                
                pos = str(i + int(originalPos))
                region = ".".join([chrom, pos])
                
                if region in genotype:
                    genotype[region].add(output)
                else:
                    genotype[region] = set([output])
    
    return(genotype)
    
    

if __name__ == "__main__":
    main()
