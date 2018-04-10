# Bioinformatics in R

## DNA-Sequence-Statistics-1
Answers to exercises in DNA Statistics (1) 
http://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter1.html

### Questions
1. What are the last twenty nucleotides of the DEN-1 Dengue virus genome sequence?
2. What is the length in nucleotides of the genome sequence for the bacterium Mycobacterium leprae strain TN (accession NC_002677)?
3. How many of each of the four nucleotides A, C, T and G, and any other symbols, are there in the Mycobacterium leprae TN genome sequence?
4. What is the GC content of the Mycobacterium leprae TN genome sequence, when (i) all non-A/C/T/G nucleotides are included, (ii) non-A/C/T/G nucleotides are discarded?
5. How many of each of the four nucleotides A, C, T and G are there in the complement of the Mycobacterium leprae TN genome sequence?
6. How many occurrences of the DNA words CC, CG and GC occur in the Mycobacterium leprae TN genome sequence?
7. How many occurrences of the DNA words CC, CG and GC occur in the (i) first 1000 and (ii) last 1000 nucleotides of the Mycobacterium leprae TN genome sequence?

## DNA-Sequence-Statistics-2
Answers to exercises in DNA Statistics (2)
http://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter2.html


### Questions
1. Draw a sliding window plot of GC content in the DEN-1 Dengue virus genome, using a window size of 200 nucleotides. Do you see any regions of unusual DNA content in the genome (eg. a high peak or low trough)?
2. Draw a sliding window plot of GC content in the genome sequence for the bacterium Mycobacterium leprae strain TN (accession NC_002677) using a window size of 20000 nucleotides. Do you see any regions of unusual DNA content in the genome (eg. a high peak or low trough)?
3. Write a function to calculate the AT content of a DNA sequence (ie. the fraction of the nucleotides in the sequence that are As or Ts). What is the AT content of the Mycobacterium leprae TN genome?
4. Write a function to draw a sliding window plot of AT content. Use it to make a sliding window plot of AT content along the Mycobacterium leprae TN genome, using a windowsize of 20000 nucleotides. Do you notice any relationship between the sliding window plot of GC content along the Mycobacterium leprae genome, and the sliding window plot of AT content?
5.  Is the 3-nucleotide word GAC GC over-represented or under-represented in the Mycobacterium leprae TN genome sequence?

## Sequence Databases
Answers to exercises in Sequence Databases 
http://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter3.html

### Questions
1. What information about the rabies virus sequence (NCBI accession NC_001542) can you obtain from its annotations in the NCBI Sequence Database?
2. How many nucleotide sequences are there from the bacterium Chlamydia trachomatis in the NCBI Sequence Database?
3. How many nucleotide sequences are there from the bacterium Chlamydia trachomatis in the RefSeq part of the NCBI Sequence Database?
4. How many nucleotide sequences were submitted to NCBI by Matthew Berriman?
5. How many nucleotide sequences from the nematode worms are there in the RefSeq Database?
6. How many nucleotide sequences for collagen genes from nematode worms are there in the NCBI Database?
7. How many mRNA sequences for collagen genes from nematode worms are there in the NCBI Database?

## Pairwise sequence aligment
Answers to exercises in sequence alignment
http://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter4.html

### Questions
1. Download FASTA-format files of the Brugia malayi Vab-3 protein (UniProt accession A8PZ80) and the Loa loa Vab-3 protein (UniProt accession E1FTG0) sequences from UniProt.
2. What is the alignment score for the optimal global alignment between the Brugia malayi Vab-3 protein and the Loa loa Vab-3 protein, when you use the BLOSUM50 scoring matrix, a gap opening penalty of -10 and a gap extension penalty of -0.5?
3. Use the printPairwiseAlignment() function to view the optimal global alignment between Brugia malayi Vab-3 protein
4. What global alignment score do you get for the two Vab-3 proteins, when you use the BLOSUM62 alignment matrix, a gap opening penalty of -10 and a gap extension penalty of -0.5?and the Loa loa Vab-3 protein, using the BLOSUM50 scoring matrix, a gap opening penalty of -10 and a gap extension penalty of -0.5.
5. What is the statistical significance of the optimal global alignment for the Brugia malayi and Loa loa Vab-3 proteins made using the BLOSUM50 scoring matrix, with a gap opening penalty of -10 and a gap extension penalty of -0.5?
6. What is the optimal global alignment score between the Brugia malayi Vab-6 protein and the Mycobacterium leprae chorismate lyase protein?

## Gene finding
Answers to exercises in gene finding
http://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter7.html#exercises

### Questions
1. How many ORFs are there on the forward strand of the DEN-1 Dengue virus genome (NCBI accession NC_001477)?
2. What are the coordinates of the rightmost (most 3’, or last) ORF in the forward strand of the DEN-1 Dengue virus genome?
3. What is the predicted protein seuence for the rightmost (most 3’, or last) ORF in the forward strand of the DEN-1 Dengue virus genome?
4. How many ORFs are there of 30 nucleotides or longer in the forward strand of the DEN-1 Dengue virus genome seuence?
5. How many ORFs longer than 248 nucleotides are there in the forward strand of the DEN-1 Dengue genome seuence?
7. How many ORFs are there on the forward strand of the rabies virus genome (NCBI accession NC_001542)?
8. What is the length of the longest ORF among the 99% of longest ORFs in 10 random seuences of the same lengths and composition as the rabies virus genome seuence?
9. How many ORFs are there in the rabies virus genome that are longer than the threshold length that you found in 8?


