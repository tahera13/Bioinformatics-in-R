#Pairwise alignment

#Set the working directory to a subfolder within the current working directory
setwd("~/Desktop/Working Directory")
setwd(paste0(getwd(), "/Statistical Analysis"))

#install biostrings
source("http://bioconductor.org/biocLite.R")
biocLite()

#libraries
library(seqinr)
library(Biostrings)

#Download FASTA-format files of the Brugia malayi Vab-3 protein (UniProt accession A8PZ80) and the Loa loa Vab-3 protein (UniProt accession E1FTG0) sequences from UniProt.
  #choose sub-database
  choosebank("swissprot")
  #fetch sequence
  bmv <- query("brugia", "AC=A8PZ80")
  closebank()
  bmvseq <- getSequence(bmv$req[[1]])
  llv <- query("loa", "AC=E1FTG0")
  llvseq <- getSequence(llv$req[[1]])
  

#What is the alignment score for the optimal global alignment between the Brugia malayi Vab-3 protein and the Loa loa Vab-3 protein, when you use the BLOSUM50 scoring matrix, a gap opening penalty of -10 and a gap extension penalty of -0.5?
  #load the coring matrix
  data("BLOSUM50")
  #convert sequence to string
  bmvseqstr <- c2s(bmvseq)
  llvseqstr <- c2s(llvseq)
  #convert str to uppercase
  bmvseqstr <- toupper(bmvseqstr)
  llvseqstr <- toupper(llvseqstr)
  #align sequences
  globalalign <- pairwiseAlignment(bmvseqstr, llvseqstr, substitutionMatrix = "BLOSUM50", gapOpening = -10, gapExtension = -0.5, scoreOnly = FALSE)
  print(globalalign)
  
#Use the printPairwiseAlignment() function to view the optimal global alignment between Brugia malayi Vab-3 protein and the Loa loa Vab-3 protein, using the BLOSUM50 scoring matrix, a gap opening penalty of -10 and a gap extension penalty of -0.5.
  #use printpairsiealignment function
  printPairwiseAlignment <- function(alignment, chunksize=60, returnlist=FALSE)
  {
    require(Biostrings)           # This function requires the Biostrings package
    seq1aln <- pattern(alignment) # Get the alignment for the first sequence
    seq2aln <- subject(alignment) # Get the alignment for the second sequence
    alnlen  <- nchar(seq1aln)     # Find the number of columns in the alignment
    starts  <- seq(1, alnlen, by=chunksize)
    n       <- length(starts)
    seq1alnresidues <- 0
    seq2alnresidues <- 0
    for (i in 1:n) {
      chunkseq1aln <- substring(seq1aln, starts[i], starts[i]+chunksize-1)
      chunkseq2aln <- substring(seq2aln, starts[i], starts[i]+chunksize-1)
      # Find out how many gaps there are in chunkseq1aln:
      gaps1 <- countPattern("-",chunkseq1aln) # countPattern() is from Biostrings package
      # Find out how many gaps there are in chunkseq2aln:
      gaps2 <- countPattern("-",chunkseq2aln) # countPattern() is from Biostrings package
      # Calculate how many residues of the first sequence we have printed so far in the alignment:
      seq1alnresidues <- seq1alnresidues + chunksize - gaps1
      # Calculate how many residues of the second sequence we have printed so far in the alignment:
      seq2alnresidues <- seq2alnresidues + chunksize - gaps2
      if (returnlist == 'FALSE')
      {
        print(paste(chunkseq1aln,seq1alnresidues))
        print(paste(chunkseq2aln,seq2alnresidues))
        print(paste(' '))
      }
    }
    if (returnlist == 'TRUE')
    {
      vector1 <- s2c(substring(seq1aln, 1, nchar(seq1aln)))
      vector2 <- s2c(substring(seq2aln, 1, nchar(seq2aln)))
      mylist <- list(vector1, vector2)
      return(mylist)
    }
  }
  #printpairwise alignment
  printPairwiseAlignment(globalalign)
  
#What global alignment score do you get for the two Vab-3 proteins, when you use the BLOSUM62 alignment matrix, a gap opening penalty of -10 and a gap extension penalty of -0.5?
  #load scoring matrix
  data("BLOSUM62")
  #align sequences
  globalalign2 <- PairwiseAlignments(bmvseqstr, llvseqstr, substitutionMatrix = "BLOSUM62", gapOpening = -10, gapExtension = -0.5, scoreOnly = FALSE)
  print(globalalign2)
  
#What is the statistical significance of the optimal global alignment for the Brugia malayi and Loa loa Vab-3 proteins made using the BLOSUM50 scoring matrix, with a gap opening penalty of -10 and a gap extension penalty of -0.5?
  #function to create a multinomial model in which the probabilities of the amino acids are set equal to their frequencies of the inputseq
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
  #use function to create random sequences with same aa frequency as inoutseq
  randomseq <- generateSeqsWithMultinomialModel(bmvseqstr, 1000)
    #align inputseq with 1000 randomn sequence
    #create numeric vector with 1000 elements
    randomscores <- double(1000)
    #run alignment function
    for (i in 1000) {
     score <- pairwiseAlignment(bmvseqstr, randomseq1[i], substitutionMatrix = "BLOSUM50", gapOpening = -9.5, gapExtension = -0.5, scoreOnly = FALSE) 
     randomscores[i] <- score
    }
  #sum of randomn scores and check if it greater than 777.5 (the score for aligning the Brugia and Loa proteins using BLOSUM50)
  print(sum(randomnscores >= 777.5))
  
#What is the optimal global alignment score between the Brugia malayi Vab-6 protein and the Mycobacterium leprae chorismate lyase protein?
  #retrieve Mycobac protein sequence
  choosebank("swissprot")
  mbp <- query("mycobacprotein", "AC=Q9CD83")
  closebank()
  mbpseq <- getSequence(mbp[[1]])
  mbpseqstr <- c2s(mbpseq)
  mbpseqstr <- toupper(mbpseqstr)
  #align sequences
  globalalign3 <- pairwiseAlignment(mbpseqstr, bmvseqstr, substitutionMatrix = "BLOSUM50", gapOpening = -9.5, gapExtension = -0.5, scoreOnly = FALSE)
  print(globalalign3)
  #calculate probability
  randomseq1 <- generateSeqsWithMultinomialModel(mbpseqstr, 1000)
  randomscores1 <- numeric(1000)
  for (i in 1000) {
    scores1 <- pairwiseAlignment(mbpseqstr, randomseq1[i], substitutionMatrix = "BLOSUM50", gapOpening = -9.5, gapExtension = -0.5, scoreOnly = FALSE)
    randomscores1[i] <- scores1
  }
  #sum of randomn scores and check if it greater than 67.5 (the score for aligning the myco and brugia proteins using BLOSUM50)
  print(sum(randomscores1 >= 67.5))