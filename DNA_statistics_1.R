#DNA sequence statistics(1) exercises

#Set the working directory to a subfolder within the current working directory
setwd("~/Desktop/Working Directory")
setwd(paste0(getwd(), "/Statistical Analysis"))

#libraries
library(seqinr)
library(ape)

#What are the last twenty nucleotides of the Dengue virus genome sequence?
  #retrieve dengue virus sequence (NC_001477)
  NC_001477 <- read.GenBank("NC_001477", as.character = "TRUE")
  # save genbank sequence as fasta format
  write.dna(NC_001477, file = "NC_001477.fasta", format = "fasta", append = FALSE, nbcol = 6, colsep = "", colw = 10)
  # convert sequence into simple string for analysis
  DEN1 <- NC_001477[[1]]
  # display last 20 of DEN1
  last20 <- tail(DEN1, 20)
  print("last 20 nucleotides are")
  print(last20)
  
#What is the length in nucleotides of the genome sequence for the bacterium Mycobacterium leprae strain TN (accession NC_002677)
  #retrieve the sequence
  NC_002677 <- read.GenBank("NC_002677", as.character = "TRUE")
  #save sequence as fasta format
  write.dna(NC_002677, file ="NC_002677", format = "fasta", append = FALSE, nbcol = 6, colsep = "", colw = 10)
  #convert sequence to simple string
  Mycobac <- NC_002677[[1]]
  #calculate length
  Mycobaclength <- length(Mycobac)
  print("Length of mycobacterium is")
  print(Mycobaclength)
  
#How many of each of the four nucleotides A, C, T and G, and any other symbols, are there in the Mycobacterium leprae TN genome sequence?
  count <- count(Mycobac,1)
  print("Nucleotide frequency is")
  print(count)
  
#What is the GC content of the Mycobacterium leprae TN genome sequence, when (i) all non-A/C/T/G nucleotides are included, (ii) non-A/C/T/G nucleotides are discarded?
  #all A/C/T/G nucleotides are included
  GCcountatgc <- GC(Mycobac)
  print(GCcountatgc)
  #all non-A/C/T/G nucleotides are included
  GCcountnonatgc <- GC(Mycobac, exact = FALSE)
  print(GCcountnonatgc)
  
#How many of each of the four nucleotides A, C, T and G are there in the complement of the Mycobacterium leprae TN genome sequence?
  Mycobaccomp <- comp(Mycobac)
  count <- count(Mycobaccomp,1)
  print("Nucleotide frequency in complement sequence is")
  print(count)
  
#How many occurrences of the DNA words CC, CG and GC occur in the Mycobacterium leprae TN genome sequence?
  count <- count(Mycobac,2)
  print("Dinucleotide frequency is")
  print(count)
  
#How many occurrences of the DNA words CC, CG and GC occur in the (i) first 1000 and (ii) last 1000 nucleotides of the Mycobacterium leprae TN genome sequence?
  # display first 1000 of Mycobac
  first1000count <- count((head(Mycobac, 1000)),2)
  print("first 1000 dinucleotide frequencies are")
  print(first1000count)
  # display last 1000 of Mycobac
  last1000count <- count((tail(Mycobac, 1000)),2)
  print("first 1000 dinucleotide frequencies are")
  print(last1000count)