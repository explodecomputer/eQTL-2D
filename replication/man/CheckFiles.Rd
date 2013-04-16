\name{CheckFiles}
\alias{CheckFiles}
\title{Check files are present}
\usage{
  CheckFiles(plink, plinkfile, probefile, intlistfile)
}
\arguments{
  \item{plink}{Path to \code{plink} executable}

  \item{plinkfile}{Path to binary plinkfile (excluding any
  suffixes)}

  \item{probefile}{Path to file with expression probe data}

  \item{intlistfile}{Path to the \code{.RData} file that
  has all the target SNPs for replication}
}
\value{
  Exits with an error if any files are missing
}
\description{
  Makes sure that all the files required for the
  replication are present
}

