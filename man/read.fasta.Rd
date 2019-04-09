\name{micRocounter}
\alias{micRocounter}
\title{
Finds microsatellite content in a fasta file
}
\usage{
read.fasta(file, minrepeats, squishy)
}
\arguments{
  \item{file}{the path to the fasta file for analysis}
  \item{minrepeats}{the minimum number of repeats required to classify as microsatellite, vector of length 5}
  \item{squishy}{a vector of length 5 describing the number of permissible "mistakes" that can be found in a microsatellite before that microsatellite is discontinued}
}
\details{
Returns a dataframe with microsatellite content grouped by monomer length
}
\examples{
read.fasta("../inst/HB.test1.fa.txt", c(6,4,3,3,3), c(1,1,1,1,1))
}
\keyword{microsatellite quantification}
