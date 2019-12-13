#!/usr/bin/env python3

# Appends a filter to a mutation map file
# Filter: mutation likely to be error because they are in a cluster
#
# Takes a mutation file, and returns pileup with 
# an extra column (or if the filter column already exist appends to it) 
# with a filter label indicating if the position is in a cluster with 
# other variatns
#
# A mutation is in a cluster if there are more or equal x mutations 
# within n bps of it

# Params: 
# 1. path to pileup file with the first two columns being 
# 3. cluster size in bps (n)
# 3. number of mutations within cluster to consider it an actual cluster (x)
# 4. name of filter to append

# Usage
# python3 filterAppend_clusterMutations.py mutation.txt 50 3 cluster


import sys
import filterAppend_functions as FA
from scipy.stats import binom

def main():
    
    pileupFile = sys.argv[1]
    neighDistance = float(sys.argv[2])
    nMutCutoff = int(sys.argv[3])
    FILTER_NOPASS = sys.argv[4]
    
    
    # Finding out number of columns in pileup file
    nCols, hasFilter = FA.getColsFilter(pileupFile)
    
    # Print pileup line for which distance is below threshold
    stream = open(pileupFile, "r")
    cache = {}
    cluster = {}
    pointerLine = 0
    highestInCluster = -1
    lastReadLine = -1
    
    # Read first position
    line = stream.readline()
    if line == "":
        sys.exit(0)
    
    line = line.rstrip().split("\t")
    lastReadLine += 1 
    highestInCluster += 1 
    cluster[lastReadLine] = line
    cache[lastReadLine] = line
    
    while True:
        
        ##
        # Reads new line or finishes
        
        # Check for new line in cache first
        if (highestInCluster + 1) in cache:
            
            line = cache[highestInCluster + 1]
            
        # Or read it from file
        else:
            
            line = stream.readline()
            if line == "":
                # Last line, dump all positions in cluster
                while pointerLine in cluster:
                    printLineWithFilter(cluster, pointerLine, nCols, nMutCutoff, hasFilter, FILTER_NOPASS)
                    pointerLine += 1
                break
            else:
                line = line.rstrip().split("\t")
                lastReadLine += 1 
                cache[lastReadLine] = line
                
        ##
        # Asses if new line in cluster or prints
        if isWithinDistance(cluster, pointerLine, line, neighDistance):
            highestInCluster = highestInCluster + 1
            cluster[highestInCluster]  = line 
        else:
            
            # Print line with filter
            printLineWithFilter(cluster, pointerLine, nCols, nMutCutoff, hasFilter, FILTER_NOPASS)
            
            # Clean cache: remove last read element
            del cache[pointerLine]
            
            # Move pointer
            pointerLine += 1
            
            # Clean cluster: remove all elements not in cluster of new pointer
            highestInCluster = purgeCluster(cluster, pointerLine, neighDistance, cache)
        
        
    
    stream.close()

def isWithinDistance(cluster, pointerLine, line, neighDistance):
    
    sameChrom = cluster[pointerLine][0] == line[0]
    withinDistance = abs(int(cluster[pointerLine][1]) - int(line[1])) <= neighDistance 
    
    return withinDistance & sameChrom


def purgeCluster(cluster, pos, neighDistance, cache):
    
    # First check if line that we'll work on is in cluster dict, if not grab it from cache
    if not pos in cluster:
        cluster[pos] = cache[pos]
    
    maxVal = -1
    for i in list(cluster.keys()):
        
        # Different chromosomes
        if cluster[pos][0] != cluster[i][0]:
            del cluster[i]
            continue
            
        if abs(int(cluster[pos][1]) - int(cluster[i][1])) > neighDistance:
            del cluster[i]
        else:
            if i > maxVal:
                maxVal = i
    
    return maxVal

def printLineWithFilter(cluster, lineToPrint, nCols, nMutCutoff, hasFilter, not_passed):
    
    linePrint = cluster[lineToPrint]
    passed = len(cluster) < nMutCutoff
    
    # Appends Filter
    if not hasFilter:
        linePrint.append("")
    linePrint = FA.getLineWithFilter(linePrint, hasFilter, passed, not_passed)
    
    # Prints 
    print(*linePrint[0:(nCols)], sep = "\t")
    
    

if __name__ == "__main__":
    main()
