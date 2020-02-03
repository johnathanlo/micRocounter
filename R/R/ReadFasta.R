library(Rcpp)
##sourceCpp("src/microsat_v4.cpp")
##source("/str_rev.R")
##source("str_comp.R")
##source("FindXmers.R")

ReadFasta2 <- function(file, xmer, minrep, tolfac){#minrepeats&squishy = vectors of values corresponding to two, three, four, five, and six-mers respectively
  if (file.exists(file) && length(minrep) == 1 && length(tolfac) == 1 && class(minrep) == "numeric" && class(tolfac) == "numeric"){
    x<-FindMS2(file, xmer, minrep, tolfac)
    return(x)
  }

  if(!file.exists(file)){
    print("File does not exist.")
    return(NULL)
  }

  if(length(minrep) != 1 || class(minrep) != "numeric"){
    print("minrepeats must be a numeric vector of length 1.")
    return(NULL)
  }

  if(length(tolfac) != 5 || class(tolfac) != "numeric"){
    print("squishy must be a numeric vector of length 1.")
    return(NULL)
  }

}

