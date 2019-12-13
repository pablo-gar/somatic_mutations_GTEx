##----------------------------------------------------
## INFO
##
## Takes one SNV pileup files, and 
## prints out the base counts of reference and alternate bases, only for those
## position that have a given read coverage and reads supporting alternate allele
##
## Rules:
##
## Author: Pablo E. Garcia-Nieto
## Date: 08/24/2017
##
## Args:
## @param[0] - path to pileup file A (a pileup  file has to have at least 8 tab-separated columns: chromosome, position, 
##                                    reference base, A counts, T counts, G counts, C counts, N counts)
## @param[1] - coverage cutoff 
## @param[2] - number of read supporting variant cutoff 
## @param[3] - MAF cutoff (lower limit)
## @param[4] - MAF cutoff (upper limit)
## @param[5] - path to ouput file
##
## @return - prints to STDOUT a tab-separated file with 5 columns:
##         - chromosome
##         - position
##         - reference base 
##         - coverage 
##         - number of reads supporting alternate
##         - number of read supporing reference
##----------------------------------------------------


import sys,os
import misc
import SNV_shared as SNV


##---------------------------------------------------
## MAIN

def main():
    
    ##--------------------
    ## Paramters read in   
    pileupFileA = sys.argv[1]
    coverageCutoffA = int(sys.argv [2])
    nCutoffA = int(sys.argv[3])
    MAFcutoffA_lower = float(sys.argv[4])
    MAFcutoffA_upper = float(sys.argv[5])
    outfile = sys.argv[6]
    ##--------------------
    
    pileupConA = open(pileupFileA, "r")
    outCon = open(outfile, "w")
    
    # Reads first lines for both files
    pileupA = SNV.getPileupInfo(pileupConA)
    
    while True:
        
        if pileupA == {}:
            break
        
        # DO stuff for similar coordinates
        SNVprint(pileupA, coverageCutoffA, nCutoffA, outCon, MAFcutoffA_upper, MAFcutoffA_lower)
        
        # Read next line for both
        pileupA = SNV.getPileupInfo(pileupConA)

    pileupConA.close()
    outCon.close()

##---------------------------------------------------


##---------------------------------------------------
## METHODS

def SNVprint(pileupA, coverageCutoffA, nCutoffA, outCon, MAFcutoffA_upper, MAFcutoffA_lower):
    
    # Gets coverage in each pileup lines and stops if they don't pass the cutoffs
    coverage = SNV.getCoverage(pileupA)
    
    if coverage < coverageCutoffA: 
        return
    
    # Gets lists of base counts sorted in decreasing order
    # sortedBaseCounts has the following structure: [Base/count list, highest first] [0 = basename, 1 = count]
    
    sortedBaseCounts = SNV.getSortedBaseCounts(pileupA)
    
    
    ref = sortedBaseCounts[0][0]
    alt = sortedBaseCounts[1][0]
    
    # Gets number of reads supporting the variant(second highest in the base/count list)
    # and stops if it does not pass the cutoffs
    refCount = sortedBaseCounts[0][1]
    n_supporting = sortedBaseCounts[1][1]
    
    passed_MAF_and_nSupporting = assertMAF_and_nSupporting(refCount, n_supporting,  nCutoffA, MAFcutoffA_upper, MAFcutoffA_lower)
    
    if not passed_MAF_and_nSupporting:
        return
    
    # Print results
    chrom = pileupA["chr"]
    pos = pileupA["pos"]
    print(chrom, pos, ref, alt, coverage, n_supporting, refCount, sep = "\t", file = outCon)
    
    return
    

def assertMAF_and_nSupporting(refCount, n_supporting, nCutoffA, MAFcutoffA_upper, MAFcutoffA_lower):
    
    ''' Asserts two filters given counts of reference and alternative allele counts,
    first that the read count for the alternative allele is greated than the given cutoff and
    second that the MAF is not greater than the given MAF cutoff.
    
    This function works with two lists of length two, the firs list has two counts
    for reference alleles and the second list has the counts for the alternative alleles.
    Positions in within each list have to correspond to each other
    
    Returns True if for one of the two counts the filters have passed'''
    
    MAF = n_supporting / (n_supporting + refCount)
    
    if n_supporting >= nCutoffA and MAF > float(MAFcutoffA_lower) and MAF < float(MAFcutoffA_upper):
        return True
    
    return False

    
    

##---------------------------------------------------


##---------------------------------------------------
## MAIN EXECUTION

if __name__ == "__main__":
    main()

##---------------------------------------------------
