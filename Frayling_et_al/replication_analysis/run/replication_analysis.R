# Joseph Powell
# Replicate epistatic signals + additional analyses for the Exeter correspondence

# Fist extract the relevant snp genotype using plink (see README)
# /path/to/plink --bfile /path/to/plink_binary_file --extract snp_list.txt --recode12 --out name_of_plink_output_data


# To run, e.g.:
# cd replication_analysis/run
# R --no-save --args /path/to/name_of_plink_output_data /path/to/probe_data.txt new_interactions_list.RData < replication_analysis.R


args        <- commandArgs(T)
plinkfile   <- args[1]
probefile   <- args[2]
intlistfile <- args[3]
outfile     <- "replication_plus.RData"

source("../R/functions.R")

CheckFiles(plinkfile, probefile, intlistfile)
geno 	<- GenoIN(plinkfile)
probes  <- ReadProbeFile(probefile)
sig     <- LoadIntList(intlistfile, plinkfile, probes)
checked <- DataChecks(probes, geno)
l       <- RunReplication(sig, checked)
predict <- GenoPrediction(sig, checked)
fit 	<- CorrectionTest(sig, checked)

newsig <- l$sig
gcm    <- l$gcm
gcs    <- l$gcs
mod    <- l$mod


save(newsig, gcm, gcs, mod, predict, fit, file=outfile)
