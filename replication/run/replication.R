# Gib Hemani
# Replicate epistatic signals

# To run, e.g.:
# cd replication/run
# R --no-save --args /path/to/binary_plink_data /path/to/probe_data.txt interaction_list.RData < replication.R

args        <- commandArgs(T)
plinkfile   <- args[1]
probefile   <- args[2]
intlistfile <- args[3]
outfile     <- "replication.RData"

source("../R/functions.R")

CheckFiles(plinkfile, probefile, intlistfile)
probes  <- ReadProbeFile(probefile)
sig     <- LoadIntList(intlistfile, plinkfile, probes)
geno    <- ExtractSNPs(sig, plinkfile)
checked <- DataChecks(probes, geno)
newsig  <- RunReplication(sig, checked)

save(newsig, file=outfile)
