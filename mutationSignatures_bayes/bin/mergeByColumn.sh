#!/bin/bash

# A list of files separated by space and merges them horizontally,
# appending columns. 
# Each files has 3 columns, for the first file all columns are kept, for the rest of the files
# only the 3rd column is appended

counter=0
while [ "$1" != "" ]
do
    
    basename=$1
    basename=${basename##*/}
    basename=${basename%.txt}
    
    if [ $counter == 0 ]
    then
        out=$(echo -e "mut\tcontext\t$basename" | cat - $1)
    else
        
        newCol=$(cut -f 3 $1)
        newCol=$(echo -e "$basename\n$newCol")
        out=$(paste <(echo "$out") <(echo "$newCol"))
        
    fi
    
    counter=counter+1
    shift
    
done

echo "$out"
