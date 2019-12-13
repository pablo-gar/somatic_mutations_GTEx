#!/usr/bin/env python3


FILTER_COL = 10 # 0 based, for pileups
FILTER_COL = 7 # 0 based, for mutation files
FILTER_PASS = "PASS"

def getColsFilter(x):
    '''
    Returns number of columns + 1 in x and true if n 
    is greater or to that number 
    '''
    
    stream = open(x, "r")
    nCols = len(stream.readline().split("\t"))
    if nCols >= FILTER_COL + 1:
        hasFilter = True
    else:
        hasFilter = False
        nCols += 1
        
    stream.close()
    
    return [nCols, hasFilter]

def getLineWithFilter(line, hasFilter, passed, filter_noPass):
    
    '''
    Return the line array with a filter added to FILTE_COL, based on the
    logicals hasFilter and passed
    
    pass and not passed labels are gotten from filter_noPass (arg) and FILTER_PASS (global)
    '''
    
    if passed and not hasFilter:
        line[FILTER_COL] = FILTER_PASS
        
    if not passed and hasFilter:
        if line[FILTER_COL] == FILTER_PASS: # overwrite a previous PASS filter
            line[FILTER_COL] = filter_noPass
        else:
            line[FILTER_COL] = ";".join([line[FILTER_COL], filter_noPass])
    
    if not passed and not hasFilter:
        line[FILTER_COL] = filter_noPass
        
    return line

