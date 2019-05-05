str_rev <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")