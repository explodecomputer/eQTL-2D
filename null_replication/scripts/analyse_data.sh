#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 10G
#SBATCH --time=01:00:00
#SBATCH --partition=test
#SBATCH --array=1-1000
#SBATCH --output=job_reports/slurm-%A_%a.out

echo "Running on ${HOSTNAME}"

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}



datafile="../data/combined"
sensnp="rs67903230"
cissnp="rs13069559"
nsim=100
minvar=0
maxvar=0.5

disc="../data/disc"
rep="../data/rep"

cischr=`grep -w $cissnp $disc.bim | cut -f 1`
echo $cischr

sd="../data/scratch/${sensnp}"



for (( sim = 1; sim <= $nsim; sim++ ))
do
	echo ${sim}

	varexp=`Rscript -e "cat(runif(1, $minvar, $maxvar))"`
	echo $varexp

	# Create phenotypes using sentinel SNP

	Rscript makephen.r ${sd}/${sensnp}_rep.raw ${sd}/${i}_${sim}_rep.fam ${varexp}
	Rscript makephen.r ${sd}/${sensnp}_disc.raw ${sd}/${i}_${sim}_disc.fam ${varexp}


	# Perform scan
	episcan -A ${sd}/${cissnp}_merge_disc.egu ${sd}/${cissnp}_disc_out -t i -F 0.00000001 -I 0.0000001 -1 1 -2 2 -f ${sd}/${i}_${sim}_disc.fam -T 4
	episcan -A ${sd}/${cissnp}_merge_rep.egu ${sd}/${cissnp}_rep_out -t i -F 0.00000001 -I 0.0000001 -1 1 -2 2 -f ${sd}/${i}_${sim}_rep.fam -T 4

	cat ${sd}/${cissnp}_disc_out[0-9]* > ${sd}/${cissnp}_disc_out
	rm ${sd}/${cissnp}_disc_out[0-9]*

	cat ${sd}/${cissnp}_rep_out[0-9]* > ${sd}/${cissnp}_rep_out
	rm ${sd}/${cissnp}_rep_out[0-9]*

	Rscript formatres.r ${sd}/${cissnp}_disc_out ${sd}/${cissnp}_rep_out ${varexp} ${i} ${sim} ${cissnp} ${sensnp} ${sd}/${i}_${sim}.rdata ${sd}/${cissnp}_exclude_disc.bim.orig

done



