#!/usr/bin/env python3
'''
FUNCTION: Remove duplicates from alignment file

NOTES:
- Files must end with .bam or .sam and must be in that format.  (The program interprets all .bam endings to be BAM files, etc)
- Ideally, we would be able to do this without the filenames even.  This will allow us to pipe it in somewhere.
- Header files must be attached to the file. 
- Files should be sorted!
'''


import random, sys, os, subprocess, argparse

__author__="Ryo Kita"
def main():

    args = run_argparse()
    rmdupFromAlignment(args.input, args.output, args.paired, args.formatIn, args.formatOut, args.samtoolsPath)

#########################
def run_argparse():
    parser = argparse.ArgumentParser(description="Remove Duplicates from Paired-End Data")
    parser.add_argument('-i', '--input', required=True,  help="Input file: It must be in BAM or SAM format")
    parser.add_argument('-o', '--output',  help="Output file: It must end in .bam or .sam.  If missing, will outputa sam file as stdout")
    parser.add_argument('-p', '--paired', dest='paired', action='store_true', help="Include -p for paired end")
    parser.add_argument('-fi', '--formatIn',  help="Format of the file.  Must be either 'bam' or 'sam'. Automatically detected if missing.")
    parser.add_argument('-fo', '--formatOut',  help="Format of the file.  Must be either 'bam' or 'sam'. Automatically detected if missing.")
    parser.add_argument('-s', '--samtoolsPath', default= os.path.join(os.path.expanduser("~"), "scripts/gtex/samtools/bin/samtools"), help="Path to samtools0.19.  Note, the newest samtools will not work.")
    args = parser.parse_args()

    if args.formatIn==None:
        ending = args.input.split(".")[-1]
        if ending not in ["bam", "sam"]:
            print("ERROR: No format detected (-fi). Specify -fi or end the input filename with '.bam' or'.sam'")
            return 0
        else: args.formatIn=ending
    
    if args.formatOut==None:
        if args.output==None or  args.output.split(".")[-1] not in ["bam", "sam"]:
            print("ERROR: No format detected (-fo). Specify -fo or end the output filename with '.bam' or '.sam'")
            return 0
        else: args.formatOut=args.output.split(".")[-1] 
    print(args)
    return args

#########################
def rmdupFromAlignment(inputFile, outputFile, paired, formatIn, formatOut, samtoolsPath):

    print("Running Remove Duplicates From Ase")
    print("OPTIONS: ")
    print("Input File ", inputFile)
    print("Output File ", outputFile)
    print("Paired End? ", paired)
    print("Input File Type: ", formatIn)
    print("Output File Type: ", formatOut)
    print("Samtools0.19 Path: ", samtoolsPath)

    #STEP 1: Randomize the quality
    tmpFileName = inputFile + ".tmp.bam"
    if formatIn=="bam":
        viewSam = subprocess.Popen((samtoolsPath, "view", "-h", inputFile), stdout=subprocess.PIPE, bufsize=1)
    else: 
        viewSam = subprocess.Popen(("cat", inputFile), stdout=subprocess.PIPE, bufsize=1)
    convertBam = subprocess.Popen((samtoolsPath, "view", "-uhS", "-o", tmpFileName, "-"), stdin=subprocess.PIPE, bufsize=1)
    for line in iter(viewSam.stdout.readline, b''):
        convertBam.stdin.write(randomize_quality(line.decode()).encode())
    viewSam.stdout.close()
    viewSam.wait()
    convertBam.stdin.close()
    convertBam.wait()

    #STEP 1.9: Does the file need to be sorted?

    #STEP 2: Remove duplicates
     #If SAM requested, then convert to sam file at the same time
    if formatOut=="bam":
        rmdup = subprocess.call((samtoolsPath, "rmdup",  tmpFileName, outputFile))
    else: 
        fileObj = open(outputFile, "w")
        rmdup = subprocess.Popen((samtoolsPath, "rmdup",  tmpFileName, "-" ), stdout = subprocess.PIPE)
        convertToSam =  subprocess.Popen((samtoolsPath, "view", "-h", "-"), stdin = rmdup.stdout, stdout = fileObj)
        rmdup.wait()
        convertToSam.wait()

    #STEP 3: Clean up
    os.system("rm " +  tmpFileName)


#########################
def randomize_quality(line, seed=False):
    '''Randomize the quality of a sam formatted line'''
    if seed:         
        #set seed to function input
        random.seed(seed)
    if line[0]!="@": #ignore, but return if header
        #convert the 4th index quality score.
        lineSplit = line.split("\t")
        lineSplit[4] = str(random.randint(0,255))
        line = "\t".join(lineSplit)
    return line

#########################
if __name__=="__main__":
    main()
#########################

