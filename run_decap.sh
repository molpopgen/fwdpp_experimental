#!/bin/bash

for i in $(seq 1 1 100)
do
    python3 msprime_decapitate.py 20000 100000 > /dev/null
    ./process_decap > /dev/null
    d=`diff -q msprime_counts.txt cpp_counts.txt`
    if [ "$d" != "" ]
    then
        echo "found a diff"
        exit 1
    fi
done    
