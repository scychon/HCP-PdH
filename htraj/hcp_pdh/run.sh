#!/bin/bash

NPD=/data/cson/PdH/htraj/hcp_pdh

cd $NPD
## declare an array variable
declare -a arrr=("r6" "r8" "r10" "r12" "r15" "r17" "r20" "r22" "r25" "r27" "r30")
declare -a arrmodel=("_h1" "_h1_surf")

# get length of an array
arrlength=${#arrr[@]}
idx=10
for (( i=0; i<${arrlength}; i++ ))
do
    mkdir $NPD/"${arrr[$i]}_h1"
    mkdir $NPD/"${arrr[$i]}_h1_surf"
    if [ -d $NPD/"${arrr[$i]}_h1" ]
    then
        cd $NPD/"${arrr[$i]}_h1"
        nohup ../../../nvtmc_eam/mc_eam -f ../param_"${arrr[$i]}".dat -r 1 > job.out 2> job.err &
    fi
    if [ -d $NPD/"${arrr[$i]}_h1_surf" ]
    then
        cd $NPD/"${arrr[$i]}_h1_surf"
        nohup ../../../nvtmc_eam/mc_eam -f ../param_"${arrr[$i]}".dat -r 1 -surf > job.out 2> job.err &
    fi
done

