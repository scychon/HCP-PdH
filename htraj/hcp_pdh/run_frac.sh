#!/bin/bash

NPD=/data/cson/PdH/htraj/hcp_pdh

cd $NPD
## declare an array variable
declare -a arrr=("r6" "r8" "r10" "r12" "r15" "r17" "r20" "r22" "r25" "r27" "r30")
declare -a arrfrac=(".2" ".4" ".6" ".8")

# get length of an array
arrlength=${#arrr[@]}
arrlenfrac=${#arrfrac[@]}
idx=0
for (( i=0; i<${arrlength}; i++ ))
do
  for (( j=0; j<${arrlenfrac}; j++ ))
  do
    mkdir $NPD/"${arrr[$i]}_h${arrfrac[$j]}"
    if [ -d $NPD/"${arrr[$i]}_h${arrfrac[$j]}" ]
    then
        cd $NPD/"${arrr[$i]}_h${arrfrac[$j]}"
        nohup ../../../nvtmc_eam/mc_eam -f ../param_"${arrr[$i]}".dat -r "${arrfrac[$j]}" > job.out 2> job.err &
    fi
    (( idx++ ))
    if (( idx==33 ))
    then
        wait
    fi
  done
done

