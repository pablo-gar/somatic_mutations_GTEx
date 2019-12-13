# Takes a table and an indexed fasta file, appends a reference sequence to
# indicated column from the indicated position at size
#
# Author: Pablo Garcia
#
# Date: 2/14/2019
#
# Requires:
#   -pysam module
#
#
# Arguments:
#   1. Reference genome file (*.fasta.gz)
#   2. Table file
#   3. column index for chromosome (0-based index)
#   4. column index for position (0-based index)
#   5. column index for where to put the reference sequence (0-based index; -1 to append at end of table)
#   6. size of reference (only odd sizes supported, the position indicated will be placed at the center)



import sys, pysam

#---------------------------------------------------#
# GLOBAL
#---------------------------------------------------#

# Total number of bases to include upstream and downstream the SNV
# it has to be an odd number (i.e. 3 == 1 upstream - SNV - 1 downstream)
CONTEXT_BASES = 5 
OFFSET = round ((CONTEXT_BASES / 2) - 0.5)


#---------------------------------------------------#
# MAIN
#---------------------------------------------------#
def main():
    
    reference, genomic_table_file, chrom_index, pos_index, append_index, CONTEXT_BASES = getArgs()
    
    OFFSET = round ((CONTEXT_BASES / 2) - 0.5)
    
    genomic_table = open(genomic_table_file, "r")
    
    for line in genomic_table:
        
        # Getting genomic info
        line = line.rstrip().split("\t")
        chrom = line[chrom_index]
        pos = int(line[pos_index])
        
        # Getting reference
        start = pos - OFFSET
        end = pos + OFFSET
        ref = openReference(reference, chrom, start, end).upper()
        
        # Appending
        if append_index == -1:
            line.append(ref)
        else:
            line[append_index] = ref
        
        print(*line, sep = "\t")
    
    genomic_table.close()

    
    
#---------------------------------------------------#
# Methdos
#---------------------------------------------------#

def getArgs():
    
    reference = sys.argv[1]
    genomic_table = sys.argv[2]
    chrom_index = int(sys.argv[3])
    pos_index = int(sys.argv[4])
    append_index = int(sys.argv[5])
    CONTEXT_BASES = int(sys.argv[6])
        
    return (reference, genomic_table, chrom_index, pos_index, append_index, CONTEXT_BASES)

def openReference(x, chrom, start, end):
    x = pysam.FastaFile(filename = x)
    region = str(chrom) + ":" + str(start) + "-" + str(end)
    return x.fetch(region = region)

if __name__ == "__main__":
    main()
