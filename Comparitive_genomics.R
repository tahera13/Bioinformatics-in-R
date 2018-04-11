#Gene finding
#could not access M ulcerans dataset, replaced it with G Gorilla and Mycobacterium leprae with Chimpanzee (ptroglodyte)

#DNA sequence statistics(2) exercises

#Set the working directory to a subfolder within the current working directory
setwd("~/Desktop/Working Directory")
setwd(paste0(getwd(), "/Statistical Analysis"))

#install biostrings
source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")
library(biomaRt)

#access and fetch data from ensembl 
gorilla <- useEnsembl(biomart = "ensembl", dataset = "ggorilla_gene_ensembl")
gorillaattr <- listAttributes(gorilla)
gorillagenes <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"), mart = gorilla)

#How many Mycobacterium ulcerans genes are there in the current version of the Ensembl Bacteria database?
  #get the gene names
  gorillagenenames <- gorillagenes[[1]]
  #calculate number of genes
  print(length(gorillagenenames))
  
#How many of the Mycobacterium ulcerans Ensembl genes are protein-coding genes?
  #create biotypes table
  gorillagenetypes <- gorillagenes[[2]]
  gorillagenetypestable <- table(gorillagenetypes)
  #get the number of protein coding genes
  gorillaprtns <- gorillagenetypestable["protein_coding"]

#How many Mycobacterium ulcerans protein-coding genes have Mycobacterium leprae orthologues?
  gorillaPT <- getBM(attributes = c("ensembl_gene_id", "ptroglodytes_homolog_ensembl_gene"), filters = "biotype", values = "protein_coding", mart = gorilla)
  gorillagenenames1 <- gorillaPT[[1]]
  gorillaPTorthologs <- gorillaPT[[2]]
  myindex <- gorillaPTorthologs != ""
  gorillagenenames2 <- gorillagenenames1[myindex]
  print(length(gorillagenenames2))
  
#How many Mycobacterium ulcerans genes have Pfam domains?(used listFilters())
  gorillapf <- getBM(attributes = c("ensembl_gene_id", "gene_biotype", "pfam"), mart = gorilla)
  gorillagenenames3 <- gorillagenespf[[1]]
  gorillapfpfamnames <- gorillagenespf[[3]]
  myindex1 <- gorillapfpfamnames != ""
  gorillagenenames4 <- gorillagenenames3[myindex1]
  print(length(gorillagenenames4))