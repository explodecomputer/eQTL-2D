#!/bin/bash

# Move all the res_2_ files to be sequential from the fist run
# There are 1959 like this
# Need to add 5380 to all of them

for (( i=$1; i<=$2; i++ ))
do

  f="results/result_2_${i}.txt.gz"

  if [ ! -f ${f} ]; then
    echo "${f} does not exist!"
    continue
  fi

  j=$(($i+5380))
  o="results/result${j}.txt.gz"

  if [ -f ${o} ]; then
    echo "${o} already exists!"
    break
  fi
  
  # echo "${f} ${o}"

  cp ${f} ${o}

done

