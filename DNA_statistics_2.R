#DNA sequence statistics(2) exercises

#Set the working directory to a subfolder within the current working directory
setwd("~/Desktop/Working Directory")
setwd(paste0(getwd(), "/Statistical Analysis"))

#libraries
library(seqinr)
library(ape)
library(rentrez)

#Draw a sliding window plot of GC content in the DEN-1 Dengue virus genome, using a window size of 200 nucleotides. Do you see any regions of unusual DNA content in the genome (eg. a high peak or low trough)?
  #retrieve dengue virus sequence (NC_001477)
  NC_001477 <- read.GenBank("NC_001477", as.character = "TRUE")
  #convert sequence into simple string for analysis
  DEN1 <- NC_001477[[1]]
  #write function for custom GC window plotting
  slidingwindowGCplot <- function(windowsize, inputseq)
  {
    GCwindow <- seq(1, length(inputseq)-windowsize, by = windowsize)
    #find length of GC window
    n <- length(GCwindow)
    #create a blank vector that is the same legth as n
    Chunks <- numeric(n)
    for (i in 1:n) {
      chunk <- inputseq[GCwindow[i]:(GCwindow[i]+windowsize-1)]
      chunkGC <- GC(chunk)
      Chunks[i] <- chunkGC
    }
    #plot graph
    plot(GCwindow, Chunks, type = "b", xlab="Nucleotide start position", ylab="GC content",main=paste("GC Plot with windowsize",windowsize))
  }

#Draw a sliding window plot of GC content in the genome sequence for the bacterium Mycobacterium leprae strain TN (accession NC_002677) using a window size of 20000 nucleotides. Do you see any regions of unusual DNA content in the genome (eg. a high peak or low trough)?
  #retrieve the sequence
  NC_002677 <- read.GenBank("NC_002677", as.character = "TRUE")
  #convert sequence to simple string
  Mycobac <- NC_002677[[1]]
  #make sliding window plot
  slidingwindowGCplot(2000, Mycobac)
  
#Write a function to calculate the AT content of a DNA sequence (ie. the fraction of the nucleotides in the sequence that are As or Ts). What is the AT content of the Mycobacterium leprae TN genome?
  ATcontent <- function(inputseq)
  {
    #get A, C, T and G content
    ATGCcontent <- count(inputseq,1)
    inputseqlength <- length(inputseq)
    Acontent <- ATGCcontent[[1]]
    Tcontent <- ATGCcontent[[2]]
    ATcontent <- (Acontent + Tcontent)/inputseqlength
    print(ATcontent)
  }

#Write a function to draw a sliding window plot of AT content. Use it to make a sliding window plot of AT content along the Mycobacterium leprae TN genome, using a windowsize of 20000 nucleotides. Do you notice any relationship between the sliding window plot of GC content along the Mycobacterium leprae genome, and the sliding window plot of AT content?
  slidingwindowATplot <- function(windowsize, inputseq)
  {
    ATwindow <- seq(1, length(inputseq)-windowsize, by = windowsize)
    n <- length(ATwindow)
    ChunkATs <- numeric(n)
    for (i in 1:n) {
      chunk1 <- inputseq[ATwindow[i]:(ATwindow[i]+windowsize-1)]
      chunkAT <- AT(chunk1)
      ChunkATs <- chunkAT[i]
    }
    plot(ATwindow, ChunkATs, type = "b", xlab="Nucleotide start position", ylab="AT content",main=paste("AT Plot with windowsize",windowsize))
  }
  slidingwindowATplot(20000, Mycobac)
  
#Is the 3-nucleotide word GAC over-represented or under-represented in the Mycobacterium leprae TN genome sequence? What is the ρ (Rho) value for this word?
  #calculate observed frequency
  tricontent <- count(Mycobac,3)
  trinucleotidesum <- sum(tricontent)
  obsfrequency <- tricontent/trinucleotidesum
  print(obsfrequency)
  #calculate expected frequency
  monocontent <- count(Mycobac,3)
  mononucleotidesum <- sum(monocontent)
  expfrequency <- ((monocontent[[1]]/mononucleotidesum) + (monocontent[[2]]/mononucleotidesum) + (monocontent[[3]]/mononucleotidesum))
  print(expfrequency)
  #use rho function in R
  rho <- rho(Mycobac, wordsize = 3)
  print(rho)
  