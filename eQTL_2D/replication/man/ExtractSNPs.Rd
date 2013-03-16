\name{ExtractSNPs}
\alias{ExtractSNPs}
\title{Read in SNPs required for replication}
\usage{
  ExtractSNPs(sig, plinkfile)
}
\arguments{
  \item{sig}{Output from \link{LoadIntList}}

  \item{plink}{Path to \code{plink} executable}

  \item{plinkfile}{Path to binary plinkfile (excluding any
  suffixes)}
}
\value{
  Returns \code{matrix} of genotype data
}
\description{
  Extracts all SNPs present in interaction list and reads
  in as 0/1/2 format matrix
}

