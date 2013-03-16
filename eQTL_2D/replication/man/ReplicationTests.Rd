\name{ReplicationTests}
\alias{ReplicationTests}
\title{Replication statistical tests}
\usage{
  ReplicationTests(geno, probes, sig, i)
}
\arguments{
  \item{geno}{\code{matrix} of genotype data}

  \item{probes}{\code{data.frame} of probes}

  \item{sig}{Output from \link{LoadIntList}}

  \item{i}{Which row of \code{sig} to run the analysis on}
}
\value{
  Returns \code{data.frame} with row \code{i} complete
}
\description{
  Calculates allele frequencies, correlation between SNPs,
  class sizes
}

