#!/bin/bash

for i in $(seq 1 1 1000)
do
    python3 msprime_decapitate.py 10000 20000 > /dev/null
    ./process_decap > /dev/null
    diff -q msprime_counts.txt cpp_counts.txt
done    
