find_xmers <- function(mon_len, x){
  valid <- c(2,3,4,5,6)
  if (mon_len == 2)
  {
    xmer <- cbind(x$Twomers$Loci, x$Twomers$Lengths, x$Twomers$Sequence, x$Twomers$Header)
  }
  if (mon_len == 3)
  {
    xmer <- cbind(x$Threemers$Loci, x$Threemers$Lengths, x$Threemers$Sequence, x$Threemers$Header)
  }
  if (mon_len == 4)
  {
    xmer <- cbind(x$Fourmers$Loci, x$Fourmers$Lengths, x$Fourmers$Sequence, x$Fourmers$Header)
  }
  if (mon_len == 5)
  {
    xmer <- cbind(x$Fivemers$Loci, x$Fivemers$Lengths, x$Fivemers$Sequence, x$Fivemers$Header)
  }
  if (mon_len == 6)
  {
    xmer <- cbind(x$Sixmers$Loci, x$Sixmers$Lengths, x$Sixmers$Sequence, x$Sixmers$Header)
  }
  if (!(mon_len%in%valid))
  {
    print("Invalid argument, monomer lengths between 2 and 6 only.")
  }

  colnames(xmer) <- c("Loci", "Lengths", "Sequence", "Header")
  xmer <- as.data.frame(xmer)
  xmer<-xmer[order(xmer$Sequence),]

  #weed out homodimers and identical sequences (complements, reverses, and reverse complements)
  xm.uniq1<-unique(xmer$Sequence)
  homodimers <- c("aa", "gg", "tt", "cc")
  xm.uniq2 <- c()
  for (seq in xm.uniq1)
  {
    if (!(seq %in% homodimers))
    {
      xm.uniq2<-c(xm.uniq2, seq)
    }
  }
  banned <- c()
  xm.uniq <- c()
  for (seq in xm.uniq2)
  {
    if (!(seq %in% banned))
    {
      seq_comp <- str_comp(seq)
      seq_rev <- str_rev(seq)
      seq_rev_comp <- str_comp(str_rev(seq))
      banned <- c(banned, seq_comp, seq_rev, seq_rev_comp)
      xm.uniq <- c(xm.uniq, seq)
    }
  }
  #reorganize and sort
  #add first column of unique xmer sequences
  xm.summary <-cbind(vector(mode = "integer", length=length(xm.uniq)),
                       vector(mode = "integer", length=length(xm.uniq)),
                       vector(mode = "list", length=length(xm.uniq)),
                       vector(mode = "numeric", length=length(xm.uniq)),
                       vector(mode = "numeric", length=length(xm.uniq)),
                       vector(mode = "numeric", length=length(xm.uniq)))
  colnames(xm.summary) <- c("Total Loci", "Total Bases", "Loci/Length/Header", "Fraction of all xmers", "Fraction of all microsats", "Fraction of whole genome")
  rownames(xm.summary) <- xm.uniq
  xm.summary <- as.data.frame(xm.summary)

  #find Total Bases, Total Loci, fraction of xmers, microsats, genome
  i <- 0
  lengths <- as.numeric(levels(xmer$Lengths))[xmer$Lengths]
  locs <- as.numeric(levels(xmer$Loci))[xmer$Loci]
  headers <- as.character(levels(xmer$Header))[xmer$Header]
  xmer.sum <- 0
  new.levels.tb <- integer(length = length(xm.uniq))
  new.levels.tl <- integer(length = length(xm.uniq))
  new.levels.ll <- vector("list", length(xm.uniq))
  new.levels.pct.xmer <- numeric(length = length(xm.uniq))
  new.levels.pct.microsats <- numeric(length = length(xm.uniq))
  new.levels.pct.genome <- numeric(length = length(xm.uniq))

  for (seq in xm.uniq)#iterate through list of unique sequences
  {
    i <- i+1
    sum <- 0
    num.loci <- 0
    loc_len <- list()

    seq_comp <- str_comp(seq)
    seq_rev <- str_rev(seq)
    seq_rev_comp <- str_comp(str_rev(seq))
    seq.indices<-grep(seq, xmer$Sequence)#obtain indices for all xmer sequences that match present iteration
    seq.indices<-c(seq.indices, grep(seq_comp, xmer$Sequence))
    seq.indices<-c(seq.indices, grep(seq_rev, xmer$Sequence))
    seq.indices<-c(seq.indices, grep(seq_rev_comp, xmer$Sequence))
    seq.indices<-unique(seq.indices)
    for (index in seq.indices)#iterate through indices of matches
    {
      num.loci <- num.loci + 1
      sum <- sum + lengths[index]#sum up their lengths
      loc_len[[num.loci]] = c(locs[index], lengths[index], headers[index])#single brackets returns indexed element as a list, double brackets returns a single indexed element itself
    }
    new.levels.ll[[i]] = loc_len
    new.levels.tb[i] <- sum
    new.levels.tl[i] <- num.loci
    new.levels.pct.genome[i] <- (sum*mon_len)/(x$`Genome Size`)
    new.levels.pct.microsats[i] <- (sum*mon_len)/(x$`Total Microsat Content`)
    xmer.sum <- xmer.sum + sum
  }
  new.levels.pct.xmer = new.levels.tb/xmer.sum
  xm.summary["Total Bases"] = new.levels.tb*mon_len
  xm.summary["Total Loci"] = new.levels.tl
  xm.summary$`Loci/Length/Header` = new.levels.ll ##not sure why this works, go back and read docs
  xm.summary["Fraction of whole genome"] = new.levels.pct.genome
  xm.summary["Fraction of all microsats"] = new.levels.pct.microsats
  xm.summary["Fraction of all xmers"] = new.levels.pct.xmer
  factor(xm.summary$`Total Bases`)
  factor(xm.summary$`Total Loci`)
  factor(xm.summary$`Fraction of whole genome`)
  factor(xm.summary$`Fraction of all microsats`)
  factor(xm.summary$`Fraction of all xmers`)

  return(xm.summary)
}
