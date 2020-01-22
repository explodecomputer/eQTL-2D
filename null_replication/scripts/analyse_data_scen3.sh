#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 10G
#SBATCH --time=01:00:00
#SBATCH --partition=test
#SBATCH --array=1-400
#SBATCH --output=job_reports/slurm-%A_%a.out

echo "Running on ${HOSTNAME}"

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}



datafile="../data/combined"
disc="../data/disc"
rep="../data/rep"

code="$1"
IFS='_' read -r -a array <<< ${code}
gene="${array[0]}"
sensnp="${array[1]}"
cissnp="${array[2]}"
sim="$2"

echo $gene
echo $sensnp
echo $cissnp
echo $sim

cischr=`grep -w $cissnp $disc.bim | cut -f 1`
echo $cischr

sd="../data/scratch/${gene}_${sensnp}_${cissnp}"

varexp=`grep -w $gene ../data/finemap.txt | grep -w $sensnp | grep -w $cissnp | head -n 1 | cut -d " " -f 3`
echo $varexp

# Create phenotypes using sentinel SNP and polygenic score

# Create polygenic score
awk '{print $1,$2,rand()}' ${sd}/polygenic.txt > ${sd}/${sim}_score.txt
plink --bfile ${sd}/polygenic --keep ${disc}list.txt --score ${sd}/${sim}_score.txt --out ${sd}/${sim}_disc_score.txt
plink --bfile ${sd}/polygenic --keep ${rep}list.txt --score ${sd}/${sim}_score.txt --out ${sd}/${sim}_rep_score.txt


if (( sim % 2 ));
then
	polyvar=0
else
	polyvar=0
fi
echo "$polyvar"

Rscript makephen.r ${sd}/${sensnp}_rep.raw ${sd}/${sim}_rep.fam ${varexp} ${sd}/${sim}_rep_score.txt.profile ${polyvar}
Rscript makephen.r ${sd}/${sensnp}_disc.raw ${sd}/${sim}_disc.fam ${varexp} ${sd}/${sim}_disc_score.txt.profile ${polyvar}


# Perform scan
episcan -A ${sd}/${cissnp}_merge_disc.egu ${sd}/${sim}_disc_out -t i -F 0.00000001 -I 0.0000001 -1 1 -2 2 -f ${sd}/${sim}_disc.fam -T 1
episcan -A ${sd}/${cissnp}_merge_rep.egu ${sd}/${sim}_rep_out -t i -F 0.00000001 -I 0.0000001 -1 1 -2 2 -f ${sd}/${sim}_rep.fam -T 1

cat ${sd}/${sim}_disc_out[0-9]* > ${sd}/${sim}_disc_out
rm ${sd}/${sim}_disc_out[0-9]*

cat ${sd}/${sim}_rep_out[0-9]* > ${sd}/${sim}_rep_out
rm ${sd}/${sim}_rep_out[0-9]*

echo "analysing"
Rscript formatres.r ${sd}/${sim}_disc_out ${sd}/${sim}_rep_out ${varexp} ${code} ${sim} ${cissnp} ${sensnp} ${sd}/${sim}.scen3.rdata ${sd}/${cissnp}_exclude_disc.bim.orig
echo "done analysing"

rm ${sd}/${sim}_*


