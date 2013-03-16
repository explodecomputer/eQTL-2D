# Gib Hemani
# Replicate epistatic signals

# To run, e.g.:
# R --no-save --args /path/to/plink /path/to/binary_plink_data /path/to/probe_data.txt /path/to/interaction_list.RData < replication.R

args       <- commandArgs(T)
plink      <- args[1]
plinkfile  <- args[2]
probefile  <- args[3]
ourfile    <- args[4]
outfile    <- "replication.RData"

source("functions.R")

CheckFiles(plink, plinkfile, probefile, ourfile)
probes <- ReadProbeFile(probefile)
sig <- LoadOurfile(ourfile, plinkfile, probes)
geno <- ExtractSNPs(sig, plink, plinkfile)
checked <- DataChecks(probes, geno)
newsig <- RunReplication(sig, checked)

save(newsig, file=outfile)
