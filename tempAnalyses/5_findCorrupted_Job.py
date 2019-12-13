#!/usr/bin/env python3
# The actual submission job

import pandas, os, subprocess, sys

outTableFile = sys.argv[1]
sraId = sys.argv[2]
fastqFile1 = sys.argv[3]
fastqFile2 = sys.argv[4]

corrupted = "False"
reason = "None"


# Check for integrety by counting lines
readNumber1 =  subprocess.Popen(" ".join(["unpigz -p 8 -c", fastqFile1, "|", "wc", "-l"]), shell = True, stdout = subprocess.PIPE).communicate()[0]
readNumber2 =  subprocess.Popen(" ".join(["unpigz -p 8 -c", fastqFile2, "|", "wc", "-l"]), shell = True, stdout = subprocess.PIPE).communicate()[0]

if readNumber1 != readNumber2:
    corrupted = "True"
    reason = "Different_read_number"
    
# Remove files if corrupted
if corrupted == "True":
    os.remove(fastqFile1)
    os.remove(fastqFile2)
    print("Corrupted and removed")

outTable = open(outTableFile, "w")
outTable.write("\t".join([sraId, corrupted, reason]) + "\n")
outTable.close
