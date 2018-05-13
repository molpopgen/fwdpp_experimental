#!/bin/bash

./fwdpp_simplify 100 10 10 0.1 10 > detail_out 2> detail_err

for i in simoutput/last_edges.*.bin
do
    echo $i
    nfile=`echo $i | sed s/edges/nodes/`
    mfile=`echo $i | sed s/edges/mutations/`
    eout=`echo $i | sed s/edges/edges_msprime/ | sed s/bin/txt/`
    nout=`echo $nfile | sed s/nodes/nodes_msprime/ | sed s/bin/txt/`
    mout=`echo $mfile | sed s/mutations/mutations_msprime/ | sed s/bin/txt/`
    
    echo $eout $nout $mfile $mout
    python3 validate_code.py $nfile $i $mfile 100 $eout $nout $mout
    #python3 check_fwd_sim_data_with_msprime.py $nfile $i 100 $eout $nout $mout

done
