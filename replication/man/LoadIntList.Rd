\name{LoadIntList}
\alias{LoadIntList}
\title{Load interaction list}
\usage{
  LoadIntList(intlistfile, plinkfile, probes)
}
\arguments{
  \item{intlistfile}{Path to the \code{.RData} file that
  has all the target SNPs for replication}
}
\value{
  Returns \code{data.frame}
}
\description{
  Loads the list of interactions (\code{.RData} file), and
  returns the subset which has SNPs and probes in common
  with the replication dataset
}

