\name{RunReplication}
\alias{RunReplication}
\title{Run replication analysis}
\usage{
  RunReplication(sig, checked)
}
\arguments{
  \item{sig}{Output from \link{LoadIntList}}

  \item{checked}{Output from \link{DataChecks}}
}
\value{
  Returns \code{data.frame} with new columns for results
  from replication data
}
\description{
  Test all SNP pairs in interaction list in replication
  dataset.
}

