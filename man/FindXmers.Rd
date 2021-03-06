\name{FindXmers}
\alias{FindXmers}
\title{
Organizes output from ReadFasta()
}
\description{
This function takes raw output from the ReadFasta()) function and produces a dataframe describing microsatellite content for one type of monomer, e.g. all twomers, threemers, fourmers, etc, aggregated by the base sequence, e.g. all atgg repeats, all tct repeats, etc.
}
\usage{
FindXmers(mon_len, micro_list)
}
\arguments{
  \item{mon_len}{A numeric value between 2 and 6 describing the length of microsatellite monomer desired
}
  \item{micro_list}{the output from ReadFasta()}
}

\details{
Sorts through the dataframe produced by ReadFasta() to find all instances of a specific monomer-length microsatellite, e.g. 2mer, 3mer, and provides some accessory statistics.
}
\value{Returns a dataframe with 6 columns and variable number of rows. Columns are Total Loci, Total Bases, Loci/Length, fraction of all xmers, fraction of all microsats, and fraction of whole genome. Rows are specific nucleotide sequences. Total loci is the total number of locations on the fasta file where a microsatellite formed from a specific nucleotide sequence is found. Total bases is the total number of bases in the fasta file. Loci/Length is a list of vectors of length 2, with one for each microsatellite; the first number in the vector describes its location in the fasta file, the second number is the number of repeats. Fraction of all xmers is a comparison of the total size, by bp, of microsatellites with a specific sequence to the total size of all microsatellites of a certain monomer length. Fraction of all microsats is the same, but with comparison to all microsatellites in the entire fasta file. Fraction of whole genome is self explanatory.
}
\examples{
x <- FindXmers(mon_len = 2, micro_list = ReadFasta(system.file("example.fasta", package = "micRocounter"), c(6,4,3,3,3), c(1,1,1,1,1)))
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Johnathan Lo and Heath Blackmon
}

\keyword{ microsatellites }
\keyword{ fasta file analysis }
