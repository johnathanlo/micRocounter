\name{read.fasta}
\alias{read.fasta}
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
Returns a list with 7 elements, where elements 1-5 are sublists with 4 elements each and elements 6-7 are integers. The five sublists are named by microsatellite monomer lengths, from twomer to sixmer, and their elements are vectors containing the locations, lengths, monomer sequence, and sequence name. Each of these vectors is the same length, and their indexes are aligned to correspond to the same data points. The 2 integer elements in the returned list are genome size and total microsatellite content, both measured in base pairs.
}
\examples{
x <- read.fasta(system.file("example.fasta", package = "micRocounter"), c(6,4,3,3,3), c(1,1,1,1,1))
}
\keyword{microsatellite quantification}
