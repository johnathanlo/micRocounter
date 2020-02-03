\name{ReadFasta2}
\alias{ReadFasta2}
\title{
Finds microsatellite content in a fasta file
}
\usage{
ReadFasta(file, xmer, minrep, tolfac)
}
\arguments{
  \item{file}{the path to the fasta file for analysis}
  \item{xmer}{the monomer length to scan for, e.g. 2 for 2mers, 3 for 3mers, etc.}
  \item{minrep}{the minimum number of repeats required to classify as microsatellite}
  \item{tolfac}{the number of permissible "mistakes" that can be found in a microsatellite before that microsatellite is discontinued}
}
\details{
Returns a list with 3 elements, where element 1 is a sublist called 'Xmers' with 4 elements and elements 2 and 3 are integers called 'Genome Size (Mbp)' and 'Total Microsat Content' respectively. The four vectors contained in 'Xmers' are called 'Loci', 'Lengths', 'Sequence', and 'SequenceNames'. Each are of identical length and their indexes correspond to the same data points, such that for a returned object called 'micro_list', micro_list$Xmers$Loci[1] is the location parameter, micro_list$Xmers$Length[1] is the number of repeated units, micro_list$Xmers$Sequence[1] is the genetic sequence of a single unit, and micro_list$Xmers$SequenceNames[1] is the scaffold name, all corresponding to a single observation in the genome file.
}
\examples{
x <- ReadFasta(system.file("example.fasta", package = "micRocounter"), xmer = 2, minrep = 6, tolfac = 0)
}
\keyword{microsatellite quantification}
