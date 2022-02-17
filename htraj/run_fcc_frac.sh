#!/bin/bash
export OMP_NUM_THREADS=1
NPD=/data/cson/PdH/htraj/fcc_pdh

cd $NPD
## declare an array variable
declare -a arrr=("r6" "r8" "r10" "r12" "r15" "r17" "r20" "r22" "r25" "r27" "r30")
declare -a arrfrac=(".2" ".25" ".3" ".4" ".5" ".6" ".75" ".8")

# get length of an array
arrlength=${#arrr[@]}
arrlenfrac=${#arrfrac[@]}
idx=0
for (( j=0; j<${arrlenfrac}; j++ ))
do
  mkdir $NPD"_h${arrfrac[$j]}"
  if [ -d $NPD"_h${arrfrac[$j]}" ]
  then
    cd $NPD"_h${arrfrac[$j]}"
    cp $NPD/param_r*.dat .
    sed -i "s/fcc_pdh/fcc_pdh${arrfrac[$j]}/g" param_r*.dat
    for (( i=0; i<${arrlength}; i++ ))
    do
      mkdir $NPD"_h${arrfrac[$j]}"/"${arrr[$i]}"
    done
  fi
done

for (( i=0; i<${arrlength}; i++ ))
do
  for (( j=0; j<${arrlenfrac}; j++ ))
  do
    if [ -d $NPD"_h${arrfrac[$j]}"/"${arrr[$i]}" ]
    then
        cd $NPD"_h${arrfrac[$j]}"/"${arrr[$i]}"
        nohup ../../../nvtmc_eam/mc_eam -f ../param_"${arrr[$i]}".dat -r "${arrfrac[$j]}" > job.out 2> job.err &
    fi
    (( idx++ ))
    if (( idx==24 ))
    then
        idx=0
        wait
    fi
  done
done

