library(Rcpp)
##sourceCpp("../microsatver3g.cpp")
##source("str_rev.r")
##source("str_comp.r")
##source("find_xmers.r")

read.fasta <- function(file, minrepeats = c(6,4,3,3,3), squishy = c(0,0,0,0,0)){#minrepeats&squishy = vectors of values corresponding to two, three, four, five, and six-mers respectively
  if (file.exists(file) && length(minrepeats) == 5 && length(squishy) == 5 && class(minrepeats) == "numeric" && class(squishy) == "numeric"){
    x<-findMS(file, minrepeats, squishy)
    return(x)
  }

  if(!file.exists(file)){
    print("File does not exist.")
    return(NULL)
  }

  if(length(minrepeats) != 5 || class(minrepeats) != "numeric"){
    print("minrepeats must be a numeric vector of length 5.")
    return(NULL)
  }

  if(length(squishy) != 5 || class(squishy) != "numeric"){
    print("squishy must be a numeric vector of length 5.")
    return(NULL)
  }

}

