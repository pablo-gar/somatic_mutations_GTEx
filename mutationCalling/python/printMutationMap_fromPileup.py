#################################
# Given a pileup file, prints out mutations that passed
# MAF, n_supporting, and coverage cutoffs
# Print a mutation map with 5 columns: chrom, pos, ref, alt, context, coverage, n_support

# Pablo E. Garcia Nieto

# Params
# - pileup file
# - eliminate triallelic sites? (True/False)
# - n_support
# - coverage
# - lower maf cutoff
# - higher maf cutoff
# - path to output file

# python3 printMutationMap_fromPileup.py /scratch/users/paedugar/somaticMutationsProject/pileups/Adipose_Subcutaneous/SRR1069097.txt 4 40 0 0.5 outMap.txt

import sys

##-----------------------------------------------------
## GLOBAL
tissueI = [3,7] # columns where the nucleotides are for each bam
tissueQual = [15,19]
cntxtI = 9

# If true the dominant base is the reference genome
dominantFromGenome = True
# If true the positions where dominant is not the reference will be ignored
ignoreIfDominantIsNotReference = False

translation_table = str.maketrans("ACTGNactgn", "TGACNtgacn")
bases = ("A", "T", "G", "C")
basesDict = {"A":0, "T":1, "G":2, "C":3}

##-----------------------------------------------------

##---------------------------------------------------
## MAIN

def main():

    ## Paramters read in   
    fileCounts = sys.argv[1]
    eliminateTriallelic = bool(sys.argv[2])
    n=int(sys.argv[3])
    coverageCutoff = int(sys.argv [4])
    MAFcutoffLower = float(sys.argv[5])
    MAFcutoff = float(sys.argv[6])
    mutFile = sys.argv [7]
    
    
    fileCounts = open(fileCounts, "r")
    mutFile = open(mutFile, "w")
    
    checked = False
    
    for posLine in fileCounts:
        
        posLine = posLine.rstrip().split("\t")
        posLine[cntxtI] = posLine[cntxtI].upper()
        
        if not checked:
            contextLength = len(posLine[cntxtI])
            middleContext = round(contextLength / 2 - 0.5)
            checked = True
        
        chrom = posLine[0]
        POS = posLine[1]
        
        # Get a list only with base counts
        varTis = posLine[ tissueI[0]:tissueI[1]  ]
        
        # Find dominant variant in the tissue showing the variation (the one with most counts)
        # And make sure that there is only one
        
        
        dom = 0 # Counts for dominant base
        domI = None # Index of dominant base
        coverage = 0
        nBasesWithReads = 0
        for a in range(len(varTis)): # loops through base counts
            varTis[a] = int(varTis[a])
            
            if(varTis[a] == 0):
                continue
            
            nBasesWithReads += 1
            coverage += varTis[a]
            
            #if varTis[a] == dom: 
            #    multiple = True
            if varTis[a] > dom:
                multiple = False
                domI = a
                dom = varTis[a]
                
                
        if eliminateTriallelic and (nBasesWithReads != 2): # Eliminate triallelic sites
            continue
        #if multiple: 
        #    continue
        if coverage < coverageCutoff: 
            continue
        if contextLength != len(posLine[cntxtI]):
            continue
        
        # Define dominant base from genome or from MAF, or ignore if dominant base is not the reference base
        if dominantFromGenome:
            genomeRef = basesDict[posLine[cntxtI][middleContext]]
            domI = genomeRef
        elif ignoreIfDominantIsNotReference:
            genomeRef = basesDict[posLine[cntxtI][middleContext]]
            if domI != genomeRef:
                continue
        else:
            posLine[cntxtI] = posLine[cntxtI][:(middleContext)] + bases[domI] + posLine[cntxtI][(middleContext+1):]
            
        # Do not contitnue if the counts of reference is 0 (important only when considering reference from the genome)
        if varTis[domI] == 0:
            continue
        
        # Counting variant sites different from dominant (only in watson strand)
        for base in range(len(varTis)):
            if domI != base:
                
                n_supporting = varTis[base]
                MAF = n_supporting / coverage
                contextLength = len(posLine[cntxtI])
                
                if n_supporting >= n and MAF <= MAFcutoff and MAF > MAFcutoffLower:
                    
                    currentContext = posLine[cntxtI]
                    printMut (mutFile, chrom, POS, bases[domI], bases[base], currentContext, coverage, n_supporting)
        
    
    
    fileCounts.close()
    mutFile.close()
    
##---------------------------------------------------


##---------------------------------------------------
## METHODS


def printMut (fileName, chrom, pos, mutFrom, mutTo, context, coverage, n_supporting):
    fileName.write(chrom + "\t" + str(pos) + "\t" + str(mutFrom) + "\t" + str(mutTo) + "\t" + context + "\t" + str(coverage) + "\t" + str(n_supporting) + "\n")
    

if __name__ == "__main__":
    main()
