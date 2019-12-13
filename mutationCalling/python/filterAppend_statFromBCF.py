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
    originalPileup = sys.argv[2]
    FILTER_NOPASS = sys.argv[3]
    cutoff_vdb = float(sys.argv[4])
    cutoff_rpb = float(sys.argv[5])
    cutoff_mqb = float(sys.argv[6])
    cutoff_bqb = float(sys.argv[7])
    cutoff_mqsb = float(sys.argv[8])
    
    
    # Finding out number of columns in pileup file
    nCols, hasFilter = FA.getColsFilter(pileupFile)
    
    # Print pileup line for which distance is below threshold
    stream = open(pileupFile, "r")
    stream_original = open(originalPileup, "r")
    for line in stream :
        
        line = line.rstrip()
        line = line.split("\t")
        
        if not hasFilter:
            line.append("")
            
        chrom = line[0]
        pos = line[1]
        vdb, rpb, mqb, bqb, mqsb = get_stats_from_original_pileup(stream_original, chrom, pos)
        
        
        # Calculate probability that this is a sequencing error
        
        passed_vdb = float(vdb) >= cutoff_vdb
        passed_rpb = float(rpb) >= cutoff_rpb
        passed_mqb = float(mqb) >= cutoff_mqb
        passed_bqb = float(bqb) >= cutoff_bqb
        passed_mqsb =float(mqsb) >= cutoff_mqsb
        
        line= FA.getLineWithFilter(line, hasFilter, passed_vdb, FILTER_NOPASS + "_variant_distance_bias")
        line= FA.getLineWithFilter(line, True, passed_rpb, FILTER_NOPASS + "_read_position_bias")
        line= FA.getLineWithFilter(line, True, passed_mqb, FILTER_NOPASS + "_mapping_quality_bias")
        line= FA.getLineWithFilter(line, True, passed_bqb, FILTER_NOPASS + "_base_quality_bias")
        line= FA.getLineWithFilter(line, True, passed_mqsb, FILTER_NOPASS + "_mapping_quality_vs_strand_bias")
        
        print(*line[0:(nCols)], sep = "\t")
    
    stream.close()
    stream_original.close()

def get_stats_from_original_pileup(stream_original, chrom, pos):
    
    while True:
        original_line = stream_original.readline()
        
        if original_line == "":
            raise ValueError("Finished reading original pileup before reaching end of mutation file\nOrdering may be inconsistent between the two files")
        
        original_line = original_line.rstrip().split("\t")
        
        if chrom == original_line[0] and pos == original_line[1]:
            break
    
    results = []
    for stat in original_line[-5:]:
        if stat == ".":
            stat = 1
        results.append(stat)
    
    return results
    
    

if __name__ == "__main__":
    main()
