##----------------------------------------------------
## INFO
##
## Takes to SNV pileup files, and for the shared posisition accross the two files
## prints out the base counts of reference and alternate bases, only for those
## position that have a given read coverage and reads supporting alternate allele
##
## Rules:
## -For a position to be considered the coverage cutoff has to be met by the two files
## -For a position to be considered the supporting-variant-read cutoff has to be met by only one of the files
##
## Author: Pablo E. Garcia-Nieto
## Date: 08/14/2017
##
## Args:
## @param[0] - path to pileup file A (a pileup  file has to have at least 8 tab-separated columns: chromosome, position, 
##                                    reference base, A counts, T counts, G counts, C counts, N counts)
## @param[1] - path to pileup file B
## @param[2] - coverage cutoff for file A
## @param[3] - coverage cutoff for file B
## @param[4] - number of read supporting variant cutoff for file A
## @param[5] - number of read supporting variant cutoff for file B
## @param[6] - MAF cutoff for file A (only considering lower than)
## @param[7] - MAF cutoff for file B (only considering lower than)
## @param[8] - path to ouput file
##
## @return - prints to STDOUT a tab-separated file with 10 columns:
##         - chromosome
##         - position
##         - reference base file A
##         - alternate base file A
##         - coverage file A
##         - number of reads supporting alternate in file A
##         - reference base file B
##         - alternate base file B
##         - coverage file B
##         - number of reads supporting alternate in file B
##         - TRUE if calls in first experiment pass the threshold, FALSE otherwise 
##         - TRUE if calls in second experiment pass the threshold, FALSE otherwise
##----------------------------------------------------


import sys,os
sys.path.append(os.path.expanduser("~/scripts/pyMods"))
import misc

##---------------------------------------------------
## MAIN

def main():
    
    ##---------------------------
    ## Paramters read in   
    pileupFileA = sys.argv[1]
    pileupFileB = sys.argv[2]
    coverageCutoffA = int(sys.argv [3])
    coverageCutoffB = int(sys.argv [4])
    nCutoffA = int(sys.argv[5])
    nCutoffB = int(sys.argv[6])
    MAFcutoffA = float(sys.argv[7])
    MAFcutoffB = float(sys.argv[8])
    outfile = sys.argv[9]
    ##---------------------------
    
    pileupConA = open(pileupFileA, "r")
    pileupConB = open(pileupFileB, "r")
    outCon = open(outfile, "w")
    
    # Reads first lines for both files
    pileupA = getPileupInfo(pileupConA)
    pileupB = getPileupInfo(pileupConB)
    
    while True:
        
        if pileupA == {} or pileupB == {}:
            break
        
        # Decide if we are in same coordinate
        [readlineA, readlineB] = getSmallerCoordinates(pileupA["chr"], pileupA["pos"], pileupB["chr"], pileupB["pos"])
        
        # If not in same coordinate read line in smaller coordinate files
        if readlineA:
            pileupA = getPileupInfo(pileupConA)
            continue
        
        if readlineB:
            pileupB = getPileupInfo(pileupConB)
            continue
        
        # DO stuff for similar coordinates
        SNVcompare(pileupA, pileupB, coverageCutoffA, coverageCutoffB, nCutoffA, nCutoffB, outCon, MAFcutoffA, MAFcutoffB)
        
        # Read next line for both
        pileupA = getPileupInfo(pileupConA)
        pileupB = getPileupInfo(pileupConB)

    pileupConA.close()
    pileupConB.close()
    outCon.close()

##---------------------------------------------------


##---------------------------------------------------
## METHODS

def getLine(fileCon):
    
    line = fileCon.readline()
    line = line.rstrip().split("\t")
    
    return line
    
def getPileupInfo(fileCon, i = 3): 
    
    line = getLine(fileCon)
    
    result = getPileupInfoFromLine(line, i)
    
    return result
    

def getPileupInfoFromLine(line, i): 
    
    if line == ['']:
        return {}
    
    result = {"chr":line[0],
            "pos":int(line[1]),
            "refBase":line[2],
            "A":int(line[i+0]),
            "T":int(line[i+1]),
            "G":int(line[i+2]),
            "C":int(line[i+3]),
            "N":int(line[i+4])}
    
    
    return result



def getSmallerCoordinates(chrA, posA, chrB, posB):
    
    ''' Returns a logical vector indicating which chromosomal coordinate from a pair is smaller,
    takes as arguments the chromosomes and position corresponding to each coordinate of the pair.
    
    @param chrA - string - chromosome in first coordinate A
    @param posA - numeric - position in first coordinate A
    @param chrB - string - chromosome in second coordinate B
    @param posB - numeric - position in first coordinate B
    
    @return - list(boolean) - of length two, if first element true then coordinates in file A
    are smaller than file B, and viceversa. If both elements false, then coordiantes are equal
    '''
    
    readlineA = False
    readlineB = False
    
    if chrA < chrB:
        readlineA = True
    
    if chrA > chrB:
        readlineB = True
    
    if chrA == chrB:
        
        if posA < posB:
            readlineA = True
        
        if posA > posB:
            readlineB = True
    
    return [readlineA, readlineB]


def SNVcompare(pileupA, pileupB, coverageCutoffA, coverageCutoffB, nCutoffA, nCutoffB, outCon, MAFcutoffA, MAFcutoffB):
    
    # Gets coverage in each pileup lines and stops if they don't pass the cutoffs
    coverage = getCoveragePair(pileupA, pileupB)
    
    if coverage[0] < coverageCutoffA or coverage[1] < coverageCutoffB: 
        return
    
    # Gets lists of base counts sorted in decreasing order
    # sortedBaseCounts has the following structure: sortedBaseCounts[0 = pileupA, 1 = pileupB] [Base/count list, highest first] [0 = basename, 1 = count]
    
    sortedBaseCounts = getSortedBaseCountsPair(pileupA, pileupB)
    
    
    ref = [sortedBaseCounts[0][0][0], sortedBaseCounts[1][0][0]]
    alt = [sortedBaseCounts[0][1][0], sortedBaseCounts[1][1][0]]
    
    # Gets number of reads supporting the variant(second highest in the base/count list)
    # and stops if it does not pass the cutoffs
    refCount = [sortedBaseCounts[0][0][1], sortedBaseCounts[1][0][1]]
    n_supporting = [sortedBaseCounts[0][1][1], sortedBaseCounts[1][1][1]] # This is alt count
    
    # If alt count is 0 it means no variant in this position, then the alt allele is set to X
    for i in range(len(n_supporting)):
        if n_supporting[i] == 0: alt[i] = "X"
    
    passed_MAF_and_nSupporting = assertMAF_and_nSupporting(refCount, n_supporting,  nCutoffA, nCutoffB, MAFcutoffA, MAFcutoffB)
    
    if "TRUE" not in passed_MAF_and_nSupporting:
        return
    
    
    # Print results
    chrom = pileupA["chr"]
    pos = pileupA["pos"]
    print(chrom, pos, ref[0], alt[0], coverage[0], n_supporting[0], refCount[0], ref[1], alt[1], coverage[1], n_supporting[1], refCount[1], 
            passed_MAF_and_nSupporting[0], passed_MAF_and_nSupporting[1],
            sep = "\t", file = outCon)
    
    return
    

def getCoveragePair(A, B):
    return [getCoverage(A), getCoverage(B)]

def getCoverage(pileupDic):
    return pileupDic["A"] \
            + pileupDic["T"] \
            + pileupDic["G"] \
            + pileupDic["C"] \
            + pileupDic["N"] \


def getSortedBaseCountsPair(A,B):
    return [getSortedBaseCounts(A), getSortedBaseCounts(B)]

def getSortedBaseCounts(pileupDic):
    toSort = dict(misc.dictSubset(pileupDic, ("A", "T", "G", "C", "N")))
    pileupDicSorted = misc.sortDictByValue(toSort, reverse = True)
    return pileupDicSorted


def assertMAF_and_nSupporting(refCount, n_supporting, nCutoffA, nCutoffB, MAFcutoffA, MAFcutoffB):
    
    ''' Asserts two filters given counts of reference and alternative allele counts,
    first that the read count for the alternative allele is greated than the given cutoff and
    second that the MAF is not greater than the given MAF cutoff.
    
    This function works with two lists of length two, the firs list has two counts
    for reference alleles and the second list has the counts for the alternative alleles.
    Positions in within each list have to correspond to each other
    
    Returns True if for one of the two counts the filters have passed'''
    
    MAF = [n_supporting[i] / (n_supporting[i] + refCount[i]) for i in range(len(refCount))]
    
    asserted = ["FALSE", "FALSE"]
    
    if n_supporting[0] >= nCutoffA and MAF[0] < float(MAFcutoffA):
        asserted[0] = "TRUE"
    if n_supporting[1] >= nCutoffB and MAF[1] < float(MAFcutoffB):
        asserted[1] = "TRUE"
    
    return asserted

    
    

##---------------------------------------------------


##---------------------------------------------------
## MAIN EXECUTION

if __name__ == "__main__":
    main()

##---------------------------------------------------
