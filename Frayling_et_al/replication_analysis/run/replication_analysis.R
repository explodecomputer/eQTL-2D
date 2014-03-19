# Joseph Powell
# Replicate epistatic signals + additional analyses for the Exeter correspondence

# To run, e.g.:
# cd replication_analysis/run
# R --no-save --args /path/to/binary_plink_data /path/to/probe_data.txt input.RData < replication_analysis.R

args        <- commandArgs(T)
plinkfile   <- args[1]
probefile   <- args[2]
intlistfile <- args[3]
outfile     <- "replication_plus.RData"

source("/R/functions.R")

CheckFiles(plinkfile, probefile, intlistfile)
geno 	<- GenoIN(plinkfile)
probes  <- ReadProbeFile(probefile)
sig     <- LoadIntList(intlistfile, plinkfile, probes)
checked <- DataChecks(probes, geno)
l       <- RunReplication(sig, checked)

newsig <- l$sig
gcm    <- l$gcm
gcs    <- l$gcs
mod    <- l$mod

save(newsig, gcm, gcs, mod, file=outfile)
