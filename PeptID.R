# MHpeptide is the quasimolecular ion M+H (single charged non-fragmented Peptide)


#To test for a peptide with a lenght of 3 or 4 AAs, load the PeptID function then run it:
  
  #PeptID (389.2033)

#To test for a peptide with a lenght of 5-AAs  
#PeptID (569.4034)

#Save the results in a variabl, for example: checkIfPeptide <- PeptID(985.5036) 

#Note that 985.5036 is M+H for Adrenorphin (sequence: YGGFMRRV)

wd <- dirname(rstudioapi::getActiveDocumentContext()$path) #extracting the current path of the folder where this R file is located
setwd(wd)  # setting wd as the working directory
# The folowing are functions needed in the code and are included in the folder "Functions"
source("Functions/Combine_Matrices_Dif_Sizes.R")     #to combine two matrices of different dimensions



PeptID <- function(MHPeptide){
Water <- 18.0106
H <- 1.0073
PP <- MHPeptide - Water - H

LL <- PP-0.005
HL <- PP+0.005


library(dplyr)
#Leucine and Isoleucine are considered as the Same AA.
AA <-  c(57.0215,71.0371,87.032,97.0528,99.0684,101.0477,103.0092,113.0841,114.0429,115.0269,
         128.0586,128.095,129.0426,131.0405,137.0589,147.0684,156.1011,163.0633,186.0793)

AASymb <- c("Glycine(G)","Alanine(A)","Serine(S)","Proline(P)","Valine(V)","Threonine(T)","Cysteine(C)","Leucine(L)/Isoleucine(I)","Asparagine(N)","Aspartic acid(D)"
            ,"Glutamine(Q)","Lysine(K)","Glutamic acid(E)","Methionine(M)","Histidine(H)","Phenylalanine(F)","Arginine(R)","Tyrosine(Y)","Tryptophan(W)")
AATable <- cbind(AASymb, AA)

#Leucine/Isoleucine are assumed to appear only 3 times in total here:
AAx2 <-  c(57.0215,71.0371,87.032,97.0528,99.0684,101.0477,103.0092,113.0841,114.0429,115.0269,
           128.0586,128.095,129.0426,131.0405,137.0589,147.0684,156.1011,163.0633,186.0793
           ,57.0215,71.0371,87.032,97.0528,99.0684,101.0477,103.0092,113.0841,113.0841,114.0429,115.0269,
           128.0586,128.095,129.0426,131.0405,137.0589,147.0684,156.1011,163.0633,186.0793)


CandidateAAall <- matrix( nrow = 1, ncol = 1)


AllCombinations <- combn(AAx2, 3)   #creates all combination of the given numbers ofs
Candidate3AA <- AllCombinations[,colSums(AllCombinations) %>% between(LL,HL)]  #keep what within the range
Candidate3AA <- t(Candidate3AA )
Candidate3AA  <- unique(t(apply(Candidate3AA, 1, sort)))
CandidateAA <- Candidate3AA

if (nrow(CandidateAA) > 0){
  for(AAmz in AA){
    for(i in 1:dim(CandidateAA )[2]) {
      for(j in 1:dim(CandidateAA )[1]) {
        if(CandidateAA [j,i]==AAmz){CandidateAA [j,i] <- AASymb[match(c(AAmz),AA)]} else CandidateAA [j,i] <- CandidateAA [j,i]}}}
  CandidateAAall <- combine.mat(CandidateAAall,CandidateAA, by = "row")
}
  

AllCombinations <- combn(AAx2, 4)   #creates all combination of the given numbers of AAs
Candidate4AA   <- AllCombinations[,colSums(AllCombinations) %>% between (LL,HL)]  #keep what within the range
Candidate4AA  <- t(Candidate4AA )
CandidateAA   <- unique(t(apply(Candidate4AA, 1, sort)))

if (nrow(CandidateAA) > 0){
  for(AAmz in AA){
    for(i in 1:dim(CandidateAA )[2]) {
      for(j in 1:dim(CandidateAA )[1]) {
        if(CandidateAA [j,i]==AAmz){CandidateAA [j,i] <- AASymb[match(c(AAmz),AA)]} else CandidateAA [j,i] <- CandidateAA [j,i]}}}
  CandidateAAall <- combine.mat(CandidateAAall,CandidateAA, by = "row")
}


#To reduce the number of possible combinations,the following Amino acids:Tryptophan,Methionine,Histidine,Tyrosine,Cysteine,
#Phenylalanine, Asparagine and arginine were asssumed to only appear twice. 
#the rest is assumed to appear 3 times. Leucine/Isoleucine only 3 times in total.

AAx3 <-  c(57.0215,71.0371,87.032,97.0528,99.0684,101.0477,103.0092,113.0841,114.0429,115.0269,
           128.0586,128.095,129.0426,131.0405,137.0589,147.0684,156.1011,163.0633,186.0793
           ,57.0215,71.0371,87.032,97.0528,99.0684,101.0477,103.0092,113.0841,113.0841,114.0429
           ,115.0269,128.0586,128.095,129.0426,131.0405,137.0589,147.0684,156.1011,163.0633,186.0793
           ,57.0215,71.0371,87.032,97.0528,99.0684,101.0477,115.0269
           ,128.0586,128.095,129.0426)


AllCombinations <- combn(AAx3, 5)   #creates all combination of the given numbers of AAs
Candidate5AA   <- AllCombinations[,colSums(AllCombinations) %>% between (LL,HL)]  #keep what within the range
Candidate5AA <- t(Candidate5AA )
CandidateAA   <- unique(t(apply(Candidate5AA, 1, sort)))

if (nrow(CandidateAA) > 0){
  for(AAmz in AA){
    for(i in 1:dim(CandidateAA )[2]) {
      for(j in 1:dim(CandidateAA )[1]) {
        if(CandidateAA [j,i]==AAmz){CandidateAA [j,i] <- AASymb[match(c(AAmz),AA)]} else CandidateAA [j,i] <- CandidateAA [j,i]}}}
  CandidateAAall <- combine.mat(CandidateAAall,CandidateAA, by = "row")
}


#To reduce the number of possible combinations,the following Amino acids:Tryptophan,Methionine,Histidine,Tyrosine,Cysteine,
#,and Phenylalanine are assummed to only appear once. Asparagine and arginine were asssumed to only appear twice. 
#the rest is assumed to appear 3 times.Leucine/Isoleucine only 3 times in total.

AAx3 <-  c(57.0215,71.0371,87.032,97.0528,99.0684,101.0477,103.0092,113.0841,113.0841,114.0429,115.0269,
           128.0586,128.095,129.0426,131.0405,137.0589,147.0684,156.1011,163.0633,186.0793
           ,57.0215,71.0371,87.032,97.0528,99.0684,101.0477,113.0841,114.0429,115.0269,
           128.0586,128.095,129.0426,156.1011
           ,57.0215,71.0371,87.032,97.0528,99.0684,101.0477,
           128.0586,128.095,129.0426,156.1011)

AllCombinations <- combn(AAx3, 6)   #creates all combination of the given numbers of AA times
Candidate6AA   <- AllCombinations[,colSums(AllCombinations) %>% between (LL,HL)]  #keep what within the range
Candidate6AA <- t(Candidate6AA )
CandidateAA   <- unique(t(apply(Candidate6AA, 1, sort)))

if (nrow(CandidateAA) > 0){
  for(AAmz in AA){
    for(i in 1:dim(CandidateAA )[2]) {
      for(j in 1:dim(CandidateAA )[1]) {
        if(CandidateAA [j,i]==AAmz){CandidateAA [j,i] <- AASymb[match(c(AAmz),AA)]} else CandidateAA [j,i] <- CandidateAA [j,i]}}}
  CandidateAAall <- combine.mat(CandidateAAall,CandidateAA, by = "row")
}


AllCombinations <- combn(AAx3, 7)   #creates all combination of the given numbers of AA times
Candidate7AA   <- AllCombinations[,colSums(AllCombinations) %>% between (LL,HL)]  #keep what within the range
Candidate7AA <- t(Candidate7AA )
CandidateAA   <- unique(t(apply(Candidate7AA, 1, sort)))

if (nrow(CandidateAA) > 0){
  for(AAmz in AA){
    for(i in 1:dim(CandidateAA )[2]) {
      for(j in 1:dim(CandidateAA )[1]) {
        if(CandidateAA [j,i]==AAmz){CandidateAA [j,i] <- AASymb[match(c(AAmz),AA)]} else CandidateAA [j,i] <- CandidateAA [j,i]}}}
  CandidateAAall <- combine.mat(CandidateAAall,CandidateAA, by = "row")
}

'
#Tryptophan,Methionine,Histidine,Tyrosine,Cysteine,
#,and Phenylalanine are assummed to only appear once here. Leucine/Isoleucine only 3 times in total.
AAx2 <-  c(57.0215,71.0371,87.032,97.0528,99.0684,101.0477,103.0092,113.0841,114.0429,115.0269,
           128.0586,128.095,129.0426,131.0405,137.0589,147.0684,156.1011,163.0633,186.0793
           ,57.0215,71.0371,87.032,97.0528,99.0684,101.0477,113.0841,113.0841,114.0429,115.0269,
           128.0586,128.095,129.0426,156.1011)
'

AllCombinations <- combn(AAx3, 8)   #creates all combination of the given numbers of AA times
Candidate8AA   <- AllCombinations[,colSums(AllCombinations) %>% between (LL,HL)]  #keep what within the range
Candidate8AA <- t(Candidate8AA )
CandidateAA   <- unique(t(apply(Candidate8AA, 1, sort)))

if (nrow(CandidateAA) > 0){
  for(AAmz in AA){
    for(i in 1:dim(CandidateAA )[2]) {
      for(j in 1:dim(CandidateAA )[1]) {
        if(CandidateAA [j,i]==AAmz){CandidateAA [j,i] <- AASymb[match(c(AAmz),AA)]} else CandidateAA [j,i] <- CandidateAA [j,i]}}}
  CandidateAAall <- combine.mat(CandidateAAall,CandidateAA, by = "row")
}


CandidateAA <- CandidateAAall

if(nrow(CandidateAA)==1){paste("Feature", MHPeptide,"might not be a Peptide.")}
else{return(CandidateAA)} 

}
