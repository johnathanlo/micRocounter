# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

FindMS2 <- function(fileloc, replen, minrepeats, tolfac) {
    .Call('_micRocounter_FindMS2', PACKAGE = 'micRocounter', fileloc, replen, minrepeats, tolfac)
}

findMS <- function(fileloc, minrepeats, tolerancefactors) {
    .Call('_micRocounter_findMS', PACKAGE = 'micRocounter', fileloc, minrepeats, tolerancefactors)
}

