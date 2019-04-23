find_xmers <- function(mon_len, micro_list){
  valid <- c(2,3,4,5,6)
  if (mon_len == 2)
  {
    xmer <- cbind(micro_list$Twomers$Loci, micro_list$Twomers$Lengths, micro_list$Twomers$Sequence, micro_list$Twomers$SequenceNames)
  }
  if (mon_len == 3)
  {
    xmer <- cbind(micro_list$Threemers$Loci, micro_list$Threemers$Lengths, micro_list$Threemers$Sequence, micro_list$Threemers$SequenceNames)
  }
  if (mon_len == 4)
  {
    xmer <- cbind(micro_list$Fourmers$Loci, micro_list$Fourmers$Lengths, micro_list$Fourmers$Sequence, micro_list$Fourmers$SequenceNames)
  }
  if (mon_len == 5)
  {
    xmer <- cbind(micro_list$Fivemers$Loci, micro_list$Fivemers$Lengths, micro_list$Fivemers$Sequence, micro_list$Fivemers$SequenceNames)
  }
  if (mon_len == 6)
  {
    xmer <- cbind(micro_list$Sixmers$Loci, micro_list$Sixmers$Lengths, micro_list$Sixmers$Sequence, micro_list$Sixmers$SequenceNames)
  }
  if (!(mon_len%in%valid))
  {
    print("Invalid argument, monomer lengths between 2 and 6 only.")
  }

  colnames(xmer) <- c("Loci", "Lengths", "Sequence", "SequenceName")
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
  colnames(xm.summary) <- c("Total Loci", "Total Bases", "Location/Length/SequenceName", "Fraction of all xmers", "Fraction of all microsats", "Fraction of whole genome")
  rownames(xm.summary) <- xm.uniq
  xm.summary <- as.data.frame(xm.summary)

  #find Total Bases, Total Loci, fraction of xmers, microsats, genome
  i <- 0
  lengths <- as.numeric(levels(xmer$Lengths))[xmer$Lengths]
  locs <- as.numeric(levels(xmer$Loci))[xmer$Loci]
  sequencenames <- as.character(levels(xmer$SequenceName))[xmer$SequenceName]
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
      loc_len[[num.loci]] = c(locs[index], lengths[index], sequencenames[index])#single brackets returns indexed element as a list, double brackets returns a single indexed element itself
    }
    new.levels.ll[[i]] = loc_len
    new.levels.tb[i] <- sum
    new.levels.tl[i] <- num.loci
    new.levels.pct.genome[i] <- (sum*mon_len)/(micro_list$`Genome Size`)
    new.levels.pct.microsats[i] <- (sum*mon_len)/(micro_list$`Total Microsat Content`)
    xmer.sum <- xmer.sum + sum
  }
  new.levels.pct.xmer = new.levels.tb/xmer.sum
  xm.summary["Total Bases"] = new.levels.tb*mon_len
  xm.summary["Total Loci"] = new.levels.tl
  xm.summary$`Location/Length/SequenceName` = new.levels.ll ##not sure why this works, go back and read docs
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
