#Gene finding

#DNA sequence statistics(2) exercises

#Set the working directory to a subfolder within the current working directory
setwd("~/Desktop/Working Directory")
setwd(paste0(getwd(), "/Statistical Analysis"))

#install biostrings
source("http://bioconductor.org/biocLite.R")
biocLite(Biostrings)
library(Biostrings)
library(seqinr)

#How many ORFs are there on the forward strand of the DEN-1 Dengue virus genome (NCBI accession NC_001477)? ANSWER: 116
  #retrieve and convert sequence to simple, uppercase string
  dengue <- read.GenBank("NC_001477", as.character = "TRUE")
  den1 <- dengue[[1]]
  den1str <- c2s(den1)
  den1str <- toupper(den1str)
  #sequence <- den1str
  #find potential start and stop codons
  findPotentialStartsAndStops <- function(sequence)
  {
    # Define a vector with the sequences of potential start and stop codons
    codons            <- c("ATG", "TAA", "TAG", "TGA")
    # Find the number of occurrences of each type of potential start or stop codon
    for (i in 1:4)
    {
      codon <- codons[i]
      # Find all occurrences of codon "codon" in sequence "sequence"
      occurrences <- matchPattern(codon, sequence)
      # Find the start positions of all occurrences of "codon" in sequence "sequence"
      codonpositions <- attr((attr(occurrences,"ranges")),"start")
      # Find the total number of potential start and stop codons in sequence "sequence"
      numoccurrences <- length(codonpositions)
      if (i == 1)
      {
        # Make a copy of vector "codonpositions" called "positions"
        positions <- codonpositions
        # Make a vector "types" containing "numoccurrences" copies of "codon"
        types <- rep(codon, numoccurrences)
      }
      else
      {
        # Add the vector "codonpositions" to the end of vector "positions":
        positions   <- append(positions, codonpositions, after=length(positions))
        # Add the vector "rep(codon, numoccurrences)" to the end of vector "types":
        types       <- append(types, rep(codon, numoccurrences), after=length(types))
      }
    }
    # Sort the vectors "positions" and "types" in order of position along the input sequence:
    indices <- order(positions)
    positions <- positions[indices]
    types <- types[indices]
    # Return a list variable including vectors "positions" and "types":
    mylist <- list(positions,types)
    return(mylist)
  }
  mylist <- findPotentialStartsAndStops(den1str)
  #find orfs
  findORFsinSeq <- function(sequence)
  {
    require(Biostrings)
    # Make vectors "positions" and "types" containing information on the positions of ATGs in the sequence:
    mylist <- findPotentialStartsAndStops(sequence)
    positions <- mylist[[1]]
    types <- mylist[[2]]
    # Make vectors "orfstarts" and "orfstops" to store the predicted start and stop codons of ORFs
    orfstarts <- numeric()
    orfstops <- numeric()
    # Make a vector "orflengths" to store the lengths of the ORFs
    orflengths <- numeric()
    # Print out the positions of ORFs in the sequence:
    # Find the length of vector "positions"
    numpositions <- length(positions)
    # There must be at least one start codon and one stop codon to have an ORF.
    if (numpositions >= 2)
    {
      for (i in 1:(numpositions-1))
      {
        posi <- positions[i]
        typei <- types[i]
        found <- 0
        while (found == 0)
        {
          for (j in (i+1):numpositions)
          {
            posj  <- positions[j]
            typej <- types[j]
            posdiff <- posj - posi
            posdiffmod3 <- posdiff %% 3
            # Add in the length of the stop codon
            orflength <- posj - posi + 3
            if (typei == "ATG" && (typej == "TAA" || typej == "TAG" || typej == "TGA") && posdiffmod3 == 0)
            {
              # Check if we have already used the stop codon at posj+2 in an ORF
              numorfs <- length(orfstops)
              usedstop <- -1
              if (numorfs > 0)
              {
                for (k in 1:numorfs)
                {
                  orfstopk <- orfstops[k]
                  if (orfstopk == (posj + 2)) { usedstop <- 1 }
                }
              }
              if (usedstop == -1)
              {
                orfstarts <- append(orfstarts, posi, after=length(orfstarts))
                orfstops <- append(orfstops, posj+2, after=length(orfstops)) # Including the stop codon.
                orflengths <- append(orflengths, orflength, after=length(orflengths))
              }
              found <- 1
              break
            }
            if (j == numpositions) { found <- 1 }
          }
        }
      }
    }
    # Sort the final ORFs by start position:
    indices <- order(orfstarts)
    orfstarts <- orfstarts[indices]
    orfstops <- orfstops[indices]
    # Find the lengths of the ORFs that we have
    orflengths <- numeric()
    numorfs <- length(orfstarts)
    for (i in 1:numorfs)
    {
      orfstart <- orfstarts[i]
      orfstop <- orfstops[i]
      orflength <- orfstop - orfstart + 1
      orflengths <- append(orflengths,orflength,after=length(orflengths))
    }
    mylist <- list(orfstarts, orfstops, orflengths)
    return(mylist)
  }
  mylist <- findORFsinSeq(den1str)
  orflengths <- mylist[[3]]                   # Find the lengths of ORFs
  print(length(orflengths))
  
#What are the coordinates of the rightmost (most 3’, or last) ORF in the forward strand of the DEN-1 Dengue virus genome?
  starts <- mylist[[1]]
  stops <- mylist[[2]]
  lastOrf <- length(orflengths)
  print(starts[lastOrf])
  print(stops[lastOrf])
  
#What is the predicted protein sequence for the rightmost (most 3’, or last) ORF in the forward strand of the DEN-1 Dengue virus genome?
  startlast <- starts[lastOrf]
  stoplast <- stops[lastOrf]
  orfvector <- den1[startlast:stoplast]
  seqinr::translate(orfvector)
  
#How many ORFs are there of 30 nucleotides or longer in the forward strand of the DEN-1 Dengue virus genome sequence?
  orflengths <- mylist[[3]]
  print(summary(orflengths >= 30))
  
#How many ORFs longer than 248 nucleotides are there in the forward strand of the DEN-1 Dengue genome sequence?
  print(summary(orflengths >= 248))
  
#How many ORFs are there on the forward strand of the rabies virus genome (NCBI accession NC_001542)?
  #retrieve and convert sequence to string and uppercase
  rabies <- read.GenBank("NC_001542", as.character = "TRUE")
  rabies <- rabies[[1]]
  rabiesstr <- c2s(rabies)
  rabiesstr <- toupper(rabiesstr)
  #calculate orfs
  mylist1 <- findPotentialStartsAndStops(rabiesstr)
  mylist1 <- findORFsinSeq(rabiesstr)
  orflengths1 <- mylist1[[3]]                   # Find the lengths of ORFs
  print(length(orflengths))
  
#What is the length of the longest ORF among the 99% of longest ORFs in 10 random sequences of the same lengths and composition as the rabies virus genome sequence?
  #use function to generate multinomials
  generateSeqsWithMultinomialModel <- function(inputsequence, X)
  {
    # Change the input sequence into a vector of letters
    require("seqinr") # This function requires the SeqinR package.
    inputsequencevector <- s2c(inputsequence)
    # Find the frequencies of the letters in the input sequence "inputsequencevector":
    mylength <- length(inputsequencevector)
    mytable <- table(inputsequencevector)
    # Find the names of the letters in the sequence
    letters <- rownames(mytable)
    numletters <- length(letters)
    probabilities <- numeric() # Make a vector to store the probabilities of letters
    for (i in 1:numletters)
    {
      letter <- letters[i]
      count <- mytable[[i]]
      probabilities[i] <- count/mylength
    }
    # Make X random sequences using the multinomial model with probabilities "probabilities"
    seqs <- numeric(X)
    for (j in 1:X)
    {
      seq <- sample(letters, mylength, rep=TRUE, prob=probabilities) # Sample with replacement
      seq <- c2s(seq)
      seqs[j] <- seq
    }
    # Return the vector of random sequences
    return(seqs)
  }
  randomseqs <- generateSeqsWithMultinomialModel(rabiesstr, 10)
  randomseqorflengths <- numeric()
  for (i in 1:10) {
    randomseq <- randomseqs[i]
    mylist2 <- findPotentialStartsAndStops(randomseq)
    mylist2 <- findORFsinSeq(randomseq)
    orflengths <- mylist2[[3]]
    randomseqorflengths <- append(randomseqorflengths, orflengths, after = length(randomseqorflengths))
  }
  #calculate the longest orf in 99% of the longest orfs
  q <- quantile(randomseqorflengths, probs=c(0.99))
  print(q)
  
#How many ORFs are there in the rabies virus genome that are longer than the threshold length that you found in Q8?
  print(summary(orflengths1 > q))
  