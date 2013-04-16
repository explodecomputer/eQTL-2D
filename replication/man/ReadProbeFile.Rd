\name{ReadProbeFile}
\alias{ReadProbeFile}
\title{Read the probefile}
\usage{
  ReadProbeFile(probefile)
}
\arguments{
  \item{probefile}{Path to file with expression probe data}
}
\value{
  Returns probe data as data.frame
}
\description{
  The probefile will be a plain text file, whitespace
  separated, with rows representing IDs and columns
  representing probes. The columns should be formatted as
  follows: \describe{ \item{\code{FID}}{Family ID (first
  row must be \code{FID})} \item{\code{IID}}{Individual ID
  column. All IID entries must be unique (first row must be
  \code{IID})} \item{\code{ILMN_xxxx}, \code{ILMN_xxxx}
  etc}{The rest of the columns are all the probes available
  in the set. Headers for the remaining columns must be the
  probe IDs.}}
}

