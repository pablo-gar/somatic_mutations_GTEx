#!/usr/bin/env python3

# Takes a pilup file, and exon boundary bed file, and returns pileup positions where
# the position is greater than a specified cutoff
# cutoff is not included in output
#
# Params: 
# 1. path to pileup file with the first two columns being 
# 2. path to bed file, this should be an exon boundary file
# 3. integer cutoff
#
# Requires :
# - bedtools
# - bash

# Usage
# python3 remove_closeToExonBoundary.py pileup.txt exonBoundaries.bed 7

import sys, subprocess

def main():
    pileupFile = sys.argv[1]
    exonBoundaryFile = sys.argv[2]
    cutoff = int(sys.argv[3])
    
    
    # Finding out number of columns in pileup file
    pileupStream = open(pileupFile, "r")
    nCols = len(pileupStream.readline().split("\t"))
    pileupStream.close()
    
    
    # Print pileup line for which distance is below threshold
    process = bed_closest(pileupFile, exonBoundaryFile)
    while True:
        
        lineDistance = process.stdout.readline().decode('UTF-8').rstrip()
        
        if lineDistance == "":
            break
        else:
            lineDistance = lineDistance.split("\t")
            dist = int(lineDistance[-1])
            if dist > cutoff:
                # Removes the copy of the coordinate
                del lineDistance[1]
                print(*lineDistance[0:(nCols)], sep = "\t")
    
    pileupStream.close()
            
    
def bed_closest(pileupFile, exonBoundaryFile):
    
    '''
    Takes a pileup file and exon boundary bed file,
    returns the process to be read line by line from 
    subprocess containing the closest distance to an
    exon boundary in each line
    '''
    
    cmd= 'bedtools closest -a <(paste <(cut -f1-2 ' + pileupFile + ') <(cut -f2- ' + pileupFile+ ') | sort -k1,1 -k2,2n) -b ' + exonBoundaryFile + ' -d -t first'
    
    process = subprocess.Popen(cmd, stdout = subprocess.PIPE, shell = True, executable = "/bin/bash")
    return (process)

if __name__ == "__main__":
    main()
