#################################
# Joins the results of  SNVsUVBAMs.py # Eliminate triallelic sites
# after all jobs from cluster finish

# Pablo E. Garcia Nieto
# 2016-05-20

# Mutation matrix from(row) to (column)
#    A   T   G   C
# 
# T  -
#
# C     -

# python3 SNV_counting.py /scratch/users/paedugar/somaticMutationsProject/pileups/Adipose_Subcutaneous/SRR1069097.txt 6 40 0.0 0.5 out1.txt out2.txt out3.txt out4.txt out5.txt

import sys,os, itertools

##-----------------------------------------------------
## GLOBAL
tissueI = [3,7] # columns where the nucleotides are for each bam
tissuePos = [10,14]
tissueQual = [15,19]
cntxtI = 9

# If true the dominant base is the reference genome
dominantFromGenome = True
# If true the positions where dominant is not the reference will be ignored
ignoreIfDominantIsNotReference = False

translation_table = str.maketrans("ACTGNactgn", "TGACNtgacn")
bases = ("A", "T", "G", "C")
basesDict = {"A":0, "T":1, "G":2, "C":3}

# list of dictionaries that will store mutatio type count with context
# 00 : T>A,  01: T>T,  02: T>G,  03:T>C,  10 : C>A, 11: C>T, 12: C>G, 13:C>C
# e.g. contextMat[0][00][ATC] is the number of T>A muse 0 2 {}
# LOOK DOWN FOR THE FUNCTION THAT CREATES THE CNTXT dictionary

#contextMat = {"00":{}, "01":{}, "02":{}, "03":{}, "10":{}, "11":{}, "12":{}, "13":{}}


##-----------------------------------------------------


##----------------------------------------------------
## PARAMETERS

## Paramters read in   

fileCounts = sys.argv[1]
n=int(sys.argv[2])
coverageCutoff = int(sys.argv [3])
MAFcutoffLower = float(sys.argv[4])
MAFcutoff = float(sys.argv[5])
outCount = sys.argv [6]
outContext = sys.argv [7]
outPosition = sys.argv [8]
mutFile = sys.argv [9]
outQuals = sys.argv[10]

# Matrix that will store mutation type count
mutMat = [
        [0,0,0,0],
        [0,0,0,0]
    ]

posMat = [
        [{},{},{},{}],
        [{},{},{},{}]
    ]

qualMat = [
        [{},{},{},{}],
        [{},{},{},{}]
    ]

for mutFrom in range(2):
    for mutTo in range(4):
        for pos in range(76):
            posMat[mutFrom][mutTo][pos] = 0
        for qual in range(45):
            qualMat[mutFrom][mutTo][qual] = 0

##----------------------------------------------------


##---------------------------------------------------
## MAIN

#def main():
    
##---------------------------------------------------


##---------------------------------------------------
## METHODS

# Function that will add counts to a mutation type in a specific context
def addCon (x,mut,con):
    if con in x[mut]:
        x[mut][con] +=1
    else:
        x[mut][con] =1
    
    return x

# Function that adds value to dictionary
def addDic (dictionary, key, value):
    if key in dictionary:
        dictionary[key] += value
    else:
        dictionary[key] = value
        

def printMut (fileName, chrom, pos, mutFrom, mutTo, context, coverage, n_supporting):
    fileName.write(chrom + "\t" + str(pos) + "\t" + str(mutFrom) + "\t" + str(mutTo) + "\t" + context + "\t" + str(coverage) + "\t" + str(n_supporting) + "\n")
    
def makeCombinations(bases, lenString, soFar = []):
    
    if(lenString == 1):
        return bases
    
    
    
    

def createContextDict(contextLen): 
    '''
    dictonary of dictionaries that will store mutatio type count with context
    00 : T>A,  01: T>T,  02: T>G,  03:T>C,  10 : C>A, 11: C>T, 12: C>G, 13:C>C
    e.g. contextMat[0][00][ATC] is the number of T>A muse 0 2 {}
    '''
    
    bases = ("A", "T", "G", "C", "N")
    if contextLen % 2 == 0:
        raise ValueError("Context length has to be an even number")
    
    middleBase = round((contextLen / 2) - 0.5)
    combinations = list(itertools.product(bases, repeat = contextLen))
    contextMat = {}
    
    for fromBase in ("0", "1"):
        if fromBase == "0":
            include = "T"
        if fromBase == "1":
            include = "C"
            
        # Getting combinations that contain as the middle base the "fromBase"
        currentCombinations = []
        for i in combinations:
            if i[middleBase] == include:
                currentCombinations.append("".join(i))
                
        # Appending those combinations to the dict
        for toBase in ("0", "1", "2", "3"):
            
            if (fromBase == "0" and toBase == "1") or (fromBase == "1" and toBase == "3"):
                continue
            
            mut = fromBase + toBase
            contextMat[mut] = {}
            for i in currentCombinations:
                contextMat[mut][i] = 0
                
    return contextMat


fileCounts = open(fileCounts, "r")
mutFile = open(mutFile, "w")

c = 1

# Will check if we are reading the new (wit position count) or old format 
hasPositions = False
hasQuals = False
checked = False
for posLine in fileCounts:
    posLine = posLine.rstrip().split("\t")
    posLine[cntxtI] = posLine[cntxtI].upper()
    
    if not checked:
        checked = True
        
        # Get the context length and created the context matrix
        contextLength = len(posLine[cntxtI])
        middleContext = round(contextLength / 2 - 0.5) # Used later
        contextMat = createContextDict(contextLength)
        
        # Check if we have infor for positions frequency and quality frequencyjhk
        if len(posLine) > tissuePos[0]:
            print(posLine)
            hasPositions = True
            if len(posLine) > tissueQual[0]:
                hasQuals = True 
            
    
    chrom = posLine[0]
    POS = posLine[1]
    
    
    # Find dominant variant in the tissue showing the variation (the one with most counts)
    # And make sure that there is only one
    varTis = posLine[ tissueI[0]:tissueI[1]  ]
    
    if hasPositions:
        VarTisPos = posLine[ tissuePos[0]:tissuePos[1]  ] # Will store the counts of positions for each base
        VarTisPosFinal = [dict() for i in range(tissuePos[1] - tissuePos[0]) ]
    if hasQuals:
        qualTisPos = posLine[ tissueQual[0]:tissueQual[1]  ] # Will store the counts of positions for each base
        qualTisPosFinal = [dict() for i in range(tissueQual[1] - tissueQual[0]) ]   
    
    dom = 0
    domI = None
    coverage = 0
    nBasesWithReads = 0
    for a in range(len(varTis)): # loops through base counts
        varTis[a] = int(varTis[a])
        
        if(varTis[a] == 0):
            continue
        else:
            nBasesWithReads += 1
        
        #Store positions
        if hasPositions:
            tempPos = VarTisPos[a].split(",")
            if not tempPos[0] == "NA":
                for posCount in tempPos:
                    pos,count = posCount.split(":")
                    VarTisPosFinal[a][int(pos)] = int(count)
        #end
        
        #Store quals
        if hasQuals:
            tempPos = qualTisPos[a].split(",")
            if not tempPos[0] == "NA":
                for posCount in tempPos:
                    pos,count = posCount.split(":")
                    qualTisPosFinal[a][int(pos)] = int(count)
        #end        
        
        coverage += varTis[a]
        if varTis[a] == dom: multiple = True
        if varTis[a] > dom:
            multiple = False
            domI = a
            dom = varTis[a]
            
            
    if nBasesWithReads != 2: # Eliminate triallelic sites
        continue
    if multiple: 
        continue
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
            
            #MAF = float(varTis[base]) / (float(varTis[base]) + float(varTis[domI]))
            n_supporting = varTis[base]
            MAF = n_supporting / coverage
            contextLength = len(posLine[cntxtI])
            
            if n_supporting >= n and MAF <= MAFcutoff and MAF > MAFcutoffLower and contextLength <= 5 :
                
                #posLine[cntxtI] = posLine[cntxtI][:(middleContext - 1)] +  bases[domI] + posLine[cntxtI][middleContext:]# Context
                originalBase = base # save base to use for positional list
                
                # If "from" A, convert "from" to T, and the "to"s to the complementary base
                if domI == 0: 
                    
                    #Prints file of mutation map before doing conversions
                    currentContext = posLine[cntxtI]
                    printMut (mutFile, chrom, POS, bases[domI], bases[base], currentContext, coverage, n_supporting)
                    
                    if base in (0,2): base += 1 # A>T; G>C
                    elif base in (1,3): base -= 1 # T>A; C>G
                    
                    mutMat[0][base] += 1
                    currentContext = posLine[cntxtI][::-1].translate(translation_table)
                    contextMat = addCon(contextMat,"0"+str(base), currentContext)
                    
                    #Appends the position where the mutations occurred
                    if hasPositions:
                        for key in VarTisPosFinal[originalBase]:
                            addDic(posMat[0][base], key, VarTisPosFinal[originalBase][key])

                    #Appends the quals where the mutations occurred
                    if hasQuals:                            
                        for key in qualTisPosFinal[originalBase]:
                            addDic(qualMat[0][base], key, qualTisPosFinal[originalBase][key])
                        
                
                # If "from" T, leave as is  
                if domI == 1: 
                    
                    #Prints file of mutation map before doing conversions
                    currentContext = posLine[cntxtI] 
                    printMut (mutFile, chrom, POS, bases[domI], bases[base], currentContext, coverage, n_supporting)
                    
                    mutMat[0][base] += 1 
                    currentContext = posLine[cntxtI] 
                    contextMat = addCon( contextMat,"0"+str(base), currentContext)
                    
                    #Appends the position where the mutations occurred 
                    if hasPositions:                
                        for key in VarTisPosFinal[originalBase]:
                            addDic(posMat[0][base], key, VarTisPosFinal[originalBase][key])

                    #Appends the quals where the mutations occurred
                    if hasQuals:                            
                        for key in qualTisPosFinal[originalBase]:
                            addDic(qualMat[0][base], key, qualTisPosFinal[originalBase][key])                       
                            
                
                # If "from" G, convert "from" to C, and the "to"s to the complementary base     
                if domI == 2: 
                    
                    #Prints file of mutation map before doing conversions
                    currentContext = posLine[cntxtI] 
                    printMut (mutFile, chrom, POS, bases[domI], bases[base], currentContext, coverage, n_supporting)
                    
                    if base in (0,2): base += 1 # A>T; G>C
                    elif base in (1,3): base -= 1 # T>A; C>G
                    mutMat[1][base] += 1
                    currentContext = posLine[cntxtI][::-1].translate(translation_table)
                    contextMat = addCon( contextMat,"1"+str(base), currentContext)
                    
                    #Appends the position where the mutations occurred 
                    if hasPositions:                    
                        for key in VarTisPosFinal[originalBase]:
                            addDic(posMat[1][base], key, VarTisPosFinal[originalBase][key])                     

                    #Appends the quals where the mutations occurred
                    if hasQuals:                            
                        for key in qualTisPosFinal[originalBase]:
                            addDic(qualMat[1][base], key, qualTisPosFinal[originalBase][key])
                                            
                # If "from" C, leave as is  
                if domI == 3: 
                    
                    #Prints file of mutation map before doing conversions
                    currentContext = posLine[cntxtI] 
                    printMut (mutFile, chrom, POS, bases[domI], bases[base], currentContext, coverage, n_supporting)
                    
                    mutMat[1][base] += 1
                    currentContext = posLine[cntxtI][::-1]
                    contextMat = addCon(contextMat,"1"+str(base), currentContext)
                    
                    #Appends the position where the mutations occurred 
                    if hasPositions:                    
                        for key in VarTisPosFinal[originalBase]:
                            addDic(posMat[1][base], key, VarTisPosFinal[originalBase][key])

                    #Appends the quals where the mutations occurred
                    if hasQuals:                            
                        for key in qualTisPosFinal[originalBase]:
                            addDic(qualMat[1][base], key, qualTisPosFinal[originalBase][key])                                                               
    
    #print (posLine)
    #if c > 4: break 
    c+=1


fileCounts.close()
mutFile.close()
# Prints results, two tables, one with context and the other without it
outCount = open(outCount, "w")
for row in range(2):
    for column in range(4):
        outCount.write(str(mutMat[row][column]))
        if column == 3:
            outCount.write("\n")
        else:
            outCount.write("\t")
outCount.close()


outContext = open(outContext, "w")
for mut in sorted(contextMat.keys()):
    for con in sorted(contextMat[mut].keys()):
        outContext.write("\t".join( [mut, con, str(contextMat[mut][con]) ]))
        outContext.write("\n")

outContext.close()


outPosition = open(outPosition, "w")
if hasPositions:
    for row in range(2):
        for column in range(4):
            outPosition.write( str(row) + str(column))  
            for key in posMat[row][column]:
                outPosition.write( "\t" + str(key) + ":" + str(posMat[row][column][key]) )
            outPosition.write("\n") 

outPosition.close()


outQuals = open(outQuals, "w")
if hasQuals:
    for row in range(2):
        for column in range(4):
            outQuals.write( str(row) + str(column)) 
            for key in qualMat[row][column]:
                outQuals.write( "\t" + str(key) + ":" + str(qualMat[row][column][key]) )
            outQuals.write("\n")    

outQuals.close()                                        
