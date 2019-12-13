#!/usr/bin/env python3

# Appends a filter to a mutation map file
# Filter: mutation likely to be sequence error
#
# Takes a pileup file, and returns pileup with 
# an extra column (or if the filter column already exist appends to it) 
# with a filter label indicating if the position has a lower than cutoff
# probability of being a sequencing error
#
# uses a conservative method that assumes all bases have a
# a min sequence quality, and calculates 1 - cdf of binom prob
# where n = coverage, x = reads supporting variant, and p = sequence quality


# Params: 
# 1. path to pileup file with the first two columns being 
# 3. float cutoff probability (filter will be appended when p >= cutoff)
# 3. minimum quality found in this pileup
# 4. name of filter to append


# Usage
# python3 filterAppend_sequenceError.py pileup.txt 0.0001 30 sequencing_error


import sys
import filterAppend_functions as FA
from scipy.stats import binom

def main():
    
    pileupFile = sys.argv[1]
    cutoff = float(sys.argv[2])
    seqQuaility = float(sys.argv[3])
    seqErrorRate = 10 ** -(seqQuaility/10) # Ilumina quality scores
    FILTER_NOPASS = sys.argv[4]
    
    
    # Finding out number of columns in pileup file
    nCols, hasFilter = FA.getColsFilter(pileupFile)
    
    # Print pileup line for which distance is below threshold
    stream = open(pileupFile, "r")
    for line in stream :
        
        line = line.rstrip()
        
        if line == "":
            break
        
        line = line.split("\t")
        if not hasFilter:
            line.append("")
        
        # Calculate probability that this is a sequencing error
        coverage = int(line[5])
        alt = int(line[6])
        prob = 1 - binom.cdf(alt, coverage, seqErrorRate)
        
        
        passed = prob < cutoff
        line= FA.getLineWithFilter(line, hasFilter, passed, FILTER_NOPASS)
        
        print(*line[0:(nCols)], sep = "\t")
    
    stream.close()

if __name__ == "__main__":
    main()
