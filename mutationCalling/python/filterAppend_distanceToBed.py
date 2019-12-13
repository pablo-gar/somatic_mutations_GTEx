#!/usr/bin/env python3

# Appends a filter to a mutation map file
# Filter: distance from closest exon boundary
#
# Takes a mutation file, a bed file, and returns pileup with 
# an extra column (or if the filter column already exist appends to it) 
# with a filter label indicating if the position is greater or equal than a specified cutoff
# for distance to closest bed element


# Params: 
# 1. path to pileup file with the first two columns being 
# 2. path to bed file, this should be an exon boundary file
# 3. integer cutoff
# 4. name of filter to append

# Requires :
# - bedtools
# - bash

# For example if you want to include a filter of distance to closest exon boundary:
# python3 filterAppend_exonBoundary.py pileup.txt exonBoundaries.bed 7 close_to_exon

# For example if you want to include for RNA edit sites
# python3 filterAppend_exonBoundary.py pileup.txt RNA_editSites.bed 1 rna_edit


import sys, subprocess
import filterAppend_functions as FA

def main():
    
    pileupFile = sys.argv[1]
    exonBoundaryFile = sys.argv[2]
    cutoff = int(sys.argv[3])
    FILTER_NOPASS = sys.argv[4]
    
    
    # Finding out number of columns in pileup file
    nCols, hasFilter = FA.getColsFilter(pileupFile)
    
    
    # Print pileup line for which distance is below threshold
    process = bed_closest(pileupFile, exonBoundaryFile)
    while True:
        
        lineDistance = process.stdout.readline().decode('UTF-8').rstrip()
        
        if lineDistance == "":
            break
        
        lineDistance = lineDistance.split("\t")
        dist = int(lineDistance[-1])
        
        # Removes the copy of the coordinate and prints output
        del lineDistance[2]
        
        passed = dist >= cutoff or (dist == -1)
        lineDistance = FA.getLineWithFilter(lineDistance, hasFilter, passed, FILTER_NOPASS)
            
        print(*lineDistance[0:(nCols)], sep = "\t")
    
            
    
    

def bed_closest(pileupFile, exonBoundaryFile):
    
    '''
    Takes a pileup file and exon boundary bed file,
    returns the process to be read line by line from 
    subprocess containing the closest distance to an
    exon boundary in each line
    '''
    
    convertMutToBed = 'awk -v FS="\\t" -v OFS="\\t" \'{$2 = $2 "\\t" $2+1; print}\' ' + pileupFile
    
    cmd= 'bedtools closest -a <(' + convertMutToBed + ' | sort -k1,1 -k2,2n) -b ' + exonBoundaryFile + ' -d -t first'
    
    process = subprocess.Popen(cmd, stdout = subprocess.PIPE, shell = True, executable = "/bin/bash")
    return (process)

if __name__ == "__main__":
    main()
