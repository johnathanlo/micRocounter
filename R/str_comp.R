str_comp <- function(x){
  x<-strsplit(x, split=NULL)
  comp <- c()
  for (base in x[[1]])
  {
    if (base == 'a')
      comp <- c(comp, 't')
    if (base == 't')
      comp <- c(comp, 'a')
    if (base == 'g')
      comp <- c(comp, 'c')
    if (base == 'c')
      comp <- c(comp, 'g')
  }
  comp <- paste(comp, collapse = '')
  return(comp)
}
  