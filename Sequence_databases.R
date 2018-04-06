#Sequence databases

#Set the working directory to a subfolder within the current working directory
setwd("~/Desktop/Working Directory")
setwd(paste0(getwd(), "/Statistical Analysis"))
 
#libraries
library(seqinr)

#What information about the rabies virus sequence (NCBI accession NC_001542) can you obtain from its annotations in the NCBI Sequence Database?
  #choose sub-database
  choosebank("refseqViruses")
  #retrieve sequence
  rb <- query("rabies", "AC=NC_001542")
  #retrieve annotations
  annots <- getAnnot(rb$req[[1]])
  print(annots[1:50])
  #close sub-database
  closebank()
  
#How many nucleotide sequences are there from the bacterium Chlamydia trachomatis in the NCBI Sequence Database?
  #choose sub-database
  choosebank("genbank")  
  #retrieve search results
  ct <- query("chlamydia", "SP=Chlamydia trachomatis")
  print(ct$nelem)
  #close sub-database
  closebank()
  
#How many nucleotide sequences are there from the bacterium Chlamydia trachomatis in the RefSeq part of the NCBI Sequence Database?
  #choose sub-database
  choosebank("refseq")  
  #retrieve search results
  ct <- query("chlamydia", "SP=Chlamydia trachomatis")
  print(ct$nelem)
  #close sub-database
  closebank()
  
#How many nucleotide sequences were submitted to NCBI by Matthew Berriman?
  #choose sub-database
  choosebank("genbank")
  #retrieve search results
  mb <- query("sequences", "AU=Berriman")
  print(mb$nelem)
  #close sub-database
  closebank()

#How many nucleotide sequences from the nematode worms are there in the RefSeq Database? 
  #choose sub-database
  choosebank("refseq")  
  #retrieve search results
  nw <- query("worms", "SP=Nematoda")
  print(nw$nelem)
  #close sub-database
  closebank()
  
#How many nucleotide sequences for collagen genes from nematode worms are there in the NCBI Database?
  choosebank("genbank")  
  #retrieve search results
  cnw <- query("collagen", "SP=Nematoda AND K=collagen")
  print(cnw$nelem)
  #close sub-database
  closebank()
  
#How many mRNA sequences for collagen genes from nematode worms are there in the NCBI Database?
  choosebank("genbank")  
  #retrieve search results
  crnw <- query("collagen", "SP=Nematoda AND K=collagen AND M=mRNA")
  print(crnw$nelem)
  #close sub-database
  closebank()
