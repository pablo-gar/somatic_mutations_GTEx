#!/bin/bash

# Appens mutation maps from different individuals into a single
# matrix that can be read by dndscv for selection analyses

for sampleFile in $@
do
    sampleName=$(basename "$sampleFile")
    sampleName="${sampleName%.*}"
    nRows=$(wc -l < $sampleFile)
    
    if [ "$nRows" == "0" ] 
    then
        continue
    fi
    
    awk -v IGNORECASE=1 -v FS="\t" -v OFS="\t" -v sample="$sampleName" '{gsub(/CHR/, ""); print sample "\t" $0}' $sampleFile | cut -f 1-5
    
    #sampleRows=$(for i in $(seq 1 $nRows); do echo $sampleName; done;)
    #out=$(paste <(echo "$sampleRows") <(cat "$sampleFile"))
    #echo "$out" | perl -pe "s/CHR//;s/chr//;s/\t[ATGCNatgcn]{2,}\t.+//"
done
