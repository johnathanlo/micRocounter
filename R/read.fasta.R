library(Rcpp)
##sourceCpp("../microsatver3g.cpp")
##source("str_rev.r")
##source("str_comp.r")
##source("find_xmers.r")

read.fasta <- function(file, minrepeats, squishy){#minrepeats&squishy = vectors of values corresponding to two, three, four, five, and six-mers respectively
  x<-findMS(file, minrepeats, squishy)
  return(x)
}
