# This is a script part of the somatic mutation caller pipelie
#
# Author: Pablo Garcia
#
# Date: 8/13/17
#
# It takes a reference genome, a bam file, and a coordinate range
# It prints out the positions in the genome for which reads have 
# a two different base calls
#
# Requires:
#   -pysam module
#   -bcftools 1.6
#       - This has to be a callable version from command line
#
#
# Arguments:
#   1. Reference genome file (*.fasta.gz)
#   2. Bam file (*.bam), if index does not exit it will be created
#   3. Chromosome name - string
#   4. Start coordinate - integer
#   5. End coordinate - integer and > Start coordinate
#   6. Bases with qualities less that this value will be ignored
#   7. Should we focus only in common SNPs (Y)? or Should ignore them (N)?
#   8. Should ignore heterozygous germline sites (Y)? or do nothing about them (N)? If Y next two arguments are required


import sys, os, pysam, subprocess, re

#---------------------------------------------------#
# GLOBAL
#---------------------------------------------------#

# Total number of bases to include upstream and downstream the SNV
# it has to be an odd number (i.e. 3 == 1 upstream - SNV - 1 downstream)
CONTEXT_BASES = 5 
OFFSET = round ((CONTEXT_BASES / 2) - 0.5)

# Shall we print position counts and qual counts?
PRINT_POS = False
PRINT_QUAL = False

#---------------------------------------------------#
# MAIN
#---------------------------------------------------#
def main():
    
    reference, bam, chrom, start, end, qualCutoff, commonSNPs, checkHet, sampleId, vcfFile = getArgs()
    
    bamPile = openAndPile(bam, chrom, start, end)
    
    ref = openReference(reference, chrom, start, end)
    
    countAndPrint(bamPile, ref, chrom, int(start), int(end), qualCutoff, commonSNPs, checkHet, sampleId, vcfFile)
    
#---------------------------------------------------#
# Methdos
#---------------------------------------------------#

def getArgs():
    
    reference = sys.argv[1]
    bam = sys.argv[2]
    chrom = sys.argv[3]
    start = sys.argv[4]
    end = sys.argv[5]
    qualCutoff = int(sys.argv[6])
    
    # Reading common SNPs argument
    
    if len(sys.argv) > 7:
        commonSNPs =  sys.argv[7].upper()
    else:
        commonSNPs = "N"
    
    if commonSNPs not in ("Y", "N"): 
        raise ValueError("Arg 7 (commonSNPs) only accepts Y/N")
    
    # Reading checking for hets argument
    
    if len(sys.argv) > 8:
        checkHet = sys.argv[8].upper()
    else:
        checkHet = "N"
        
    if checkHet not in ("Y", "N"):
        raise ValueError("Arg 8 (checking for hets) only accepts Y/N")
    
    if checkHet == "Y":
        if not len(sys.argv) > 10:
            raise ValueError("Arguments 9 (sample Id) and 10 (vcfFile) are needed for checking heterozygous status")
        
        checkHet = True
        sampleId = sys.argv[9]
        vcfFile = sys.argv[10]
    else:
        checkHet = False
        sampleId = ""
        vcfFile = ""
        
        
        
    return (reference, bam, chrom, start, end, qualCutoff, commonSNPs, checkHet, sampleId, vcfFile)


def openAndPile(x, chrom, start, end):
    x = pysam.AlignmentFile (x, "rb")
    region = str(chrom) + ":" + str(start) + "-" + str(end)
    return x.pileup(chrom, region = region)

def openReference(x, chrom, start, end):
    x = pysam.FastaFile(filename = x)
    region = str(chrom) + ":" + str(start) + "-" + str(end)
    return x.fetch(region = region)
    
    
    
def countAndPrint(bamIter, ref, chrom, start, end, qualCutoff, commonSNPs, checkHet, sampleId, vcfFile):
    
    
    pileCounts, previousBase = updateCounts()
    
    
    # Goes through all positions that have a mapped read
    for pile in bamIter:
        
        refBase = ""
        isSNV = False
        
        # Goes through all read mapping in this position
        #iter = 1
        for pileRead in pile.pileups:
            
            if not refBase:
                
                # Gets refbase and context
                refPos = pile.reference_pos # 0 - based coordinates
                refBase = ref[ refPos-start - OFFSET + 1 : refPos - start + OFFSET + 2 ] # 0 based coordinates
                refBase = refBase.upper()
                
                
            #------------------------------------------------------------------------------
            # Gathering info and stoping or skiping when needed
            
            if len(refBase) != CONTEXT_BASES: break
            
            #  Igonorig or exclusively focusing in commonSNPs
            if commonSNPs == "Y":
                if refBase[OFFSET] != "N": break
            else:
                if refBase[OFFSET] == "N": break
            
            # Only work if there is an actual base mapping at this position 
            if pileRead.is_del or pileRead.is_refskip: continue
            
            queryPos = pileRead.query_position
            quality = pileRead.alignment.query_qualities[queryPos]
            if quality < qualCutoff: continue
            #------------------------------------------------------------------------------
            
            
            # Reads mapping base in read and adds it to its corresponding counter
            queryBase = pileRead.alignment.query_sequence[queryPos]
            pileCounts["counts"][queryBase] += 1
            
            # Store what position in the read the mapped base is and adds to counter that keeps track
            # of the frequency a position contains a variant
            if PRINT_POS:
                addPosition(pileCounts, pileRead.alignment, queryPos, queryBase)
            if PRINT_QUAL:
                addQuality(pileCounts, queryBase, quality)
            
            # check whehter the current base is different from the previous ones
            if not isSNV: [previousBase, isSNV] = SNVeval(previousBase, queryBase)
            
            #Turns on flag of base in group (when countin more than bam file at a time
            #if not pileGroupsCounts["flag"]: pileGroupsCounts["flag"] = 1
            
        
        # Printing out only when there is a refbase and a snv (if counting multiple files then you have to check that there is info for all of them)
        if refBase and isSNV:
            
            if checkHet:
                passedGenotype = assertGenotype(sampleId, vcfFile, chrom, refPos)
            else:
                passedGenotype = True
                
            if passedGenotype:
                
                # If there is a common SNP in the context tries to get the corresponding allele
                if "N" in refBase and vcfFile != "":
                    refBase = correctContext(refBase, sampleId, vcfFile, chrom, refPos)
                    
                # Printing results
                printResults(chrom, refPos, refBase, pileCounts, PRINT_POS, PRINT_QUAL)
                
                
        pileCounts, previousBase = updateCounts()
            
            
                                            

def correctContext(refBase, sampleId, vcfFile, chrom, refPos):
    
    '''
    Given a genomic sequence and a vcfFile, it will get assigned the correct genotype
    where possible in the Ns
    
    refPos is 0-based
    vcfFile is 1-based
    
    @return corrected sequence
    '''
    
    genotypes = getGenotypeInfo(sampleId, vcfFile, re.sub(r'chr', '', chrom) +  ":" + str(refPos - OFFSET + 1) + "-" + str(refPos + OFFSET + 1))
    
    while(len(genotypes) > 0):
        genotype = genotypes.pop()
        
        # Parses out genotype info out of bcftools output
        genotype = genotype.split("|")
        pos = int(genotype[0])
        refAlt = genotype[1:3] #Alleles
        genotype = genotype[3].split("/")
        
        if genotype[0] != genotype[1]: continue
        if genotype[0] not in ["0","1"]: continue
        
        genotype = [int(i) for i in genotype]
        allele =  refAlt[genotype[0]] # Only one allele
        
        if len(allele) > 1: continue
        
        c = pos - refPos + 1
        refBase = "".join([refBase[0:c], allele, refBase[c + 1:]])
        
    return refBase

def getGenotypeInfo(sampleId, vcfFile, region):
    
    '''
    Returns the alleles and genotype of a given sample given a region
    
    requires a callabe bcftools command
    
    @param region - format is the same as bcfTools [chr:pos][chr:start-end], ...
    '''
    
    bcfToolsCmd = ["bcftools", "query", "-f", "%POS|%REF|%ALT|[%GT]\n", "-s", sampleId, "-r", region, vcfFile]
    bcfTools = subprocess.Popen(bcfToolsCmd, stdout = subprocess.PIPE)
    genotype = bcfTools.communicate()[0].decode("ascii")
    genotype = genotype.split("\n")
    genotype = genotype[:-1]
    
    return genotype

def getGenotype(sampleId, vcfFile, region):
    
    genotype = getGenotypeInfo(sampleId, vcfFile, region)
    if len(genotype) > 0:
        genotype = [i.split("|")[-1:][0] for i in genotype]# last element
    else:
        genotype = [""]
        
    return genotype

def assertGenotype(sampleId, vcfFile, chrom, refPos):
    
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
    
    genotypes = getGenotype(sampleId, vcfFile, ":".join([re.sub(r'chr', '', chrom), str(refPos + 1)]))
    for i in genotypes:
        if i == "0/1" or i == "1/0" or i == "1/." or i == "0/." or i == "./1" or i == "./0":
            assertion = False
    
    return assertion
    
    
def updateCounts():
    a = {"counts":{"A":0,"T":0,"G":0,"C":0, "N":0, "R":0}, "flag":0, "pos":{"A":{},"T":{},"G":{},"C":{}, "N":{}}, "qual":{"A":{},"T":{},"G":{},"C":{}, "N":{}}}
    b = ""
    return [a,b]



def addPosition(counts, readAlignment, pos, queryBase):
    queryPosOriginal = getPositionInRead(readAlignment, pos)
    
    if queryPosOriginal in counts["pos"][queryBase]:
        counts["pos"][queryBase][queryPosOriginal] += 1 
    else:
        counts["pos"][queryBase][queryPosOriginal] = 1 

def addQuality(counts, base, qual):
    if qual in counts["qual"][base]:
        counts["qual"][base][qual] += 1
    else:
        counts["qual"][base][qual] = 1
    
def getPositionInRead(readAlignment, pos):
    if readAlignment.is_reverse:
        return len(readAlignment.query_sequence) - pos + 1
    else:
        return pos + 1

def SNVeval(prev, current):
    state = False
    if not current == "N":
        if not prev:
            prev = current
        else:
            if not prev == current:
                state = True
    
    return [prev, state]


def printResults(chrom, refPos, refBase, pileCounts, printPos, printQual):

    # Print refbase
    print(chrom, refPos + 1 , refBase[OFFSET], sep = "\t", end= "\t")
    pileCounts["counts"]["R"] = pileCounts["counts"][refBase[OFFSET]]
    
    # Gets string to print
    bamPileStrings = getStringToPrint(pileCounts, refBase, printPos, printQual) #refBase added just for context 

    print(bamPileStrings)

def getStringToPrint(pileCounts, refBase, printPos = True, printQual = True):
    
    # Vars to be used
    baseString = []
    posString = []
    qualString = []
    
    # Making the string with base counts
    for base in ["A", "T", "G", "C", "N", "R"]:
        baseString.append(str(pileCounts["counts"][base]))  
        if base == "R": continue
    
        # Making the string with position in the read   
        if printPos:
            posCountList = []
            if pileCounts["pos"][base]:
                for pos in pileCounts["pos"][base]:
                    posCountList.append(":".join([str(pos), str(pileCounts["pos"][base][pos])]))
                posString.append(",".join(posCountList))
            else:
                posString.append("NA")
    
        # Making the string with qualities
        if printQual:
            qualCountList = []
            if pileCounts["qual"][base]:
                for qual in pileCounts["qual"][base]:
                    qualCountList.append(":".join([str(qual), str(pileCounts["qual"][base][qual])]))
                qualString.append(",".join(qualCountList))
            else:
                qualString.append("NA")
    
    baseString = "\t".join(baseString)
    final = baseString
    
    # Append context
    final = "\t".join([final, refBase])
    
    if(printPos):
        posString = "\t".join(posString)
        final = "\t".join([final, posString])

    if(printQual):
        qualString = "\t".join(qualString)
        final = "\t".join([final, qualString])

    return final

    
#---------------------------------------------------#

if __name__ == "__main__": 
    main()
