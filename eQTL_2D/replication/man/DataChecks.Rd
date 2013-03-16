\name{DataChecks}
\alias{DataChecks}
\title{Check probe and genotype data}
\usage{
  DataChecks(probes, geno)
}
\arguments{
  \item{probes}{Output from \link{ReadProbeFile}}

  \item{geno}{Output from \link{ExtractSNPs}}
}
\value{
  Returns \code{list} of geno \code{matrix} and probes
  \code{data.frame}
}
\description{
  Checks that IIDs in geno and probes are matched, and
  reorders accordingly.
}

