#This script is for running the bsseq tool for one data set (=genomic window) at a time

#Setting up the environment 

.libPaths("/scratch/work/hallav1/R")
Sys.setenv(R_LIBS="/scratch/work/hallav1/R")


#Arguments to the script

input_args <- commandArgs(TRUE)
#The input arguments should be
#1: Folder to be set as environment
#2: Input file name (proportion table)
#3: Input file name (design matrix) NOT ACTUALLY NEEDED
#4: Output file name (+ path if different than environment)
#5: Number of columns in proportion table = total number of replicates


#Test that the number of arguments is correct
if (length(input_args)!=5){
  stop("The number of input arguments is incorrect.",call.=FALSE)
}

setwd(input_args[1])

file_id=unlist(strsplit(input_args[2],"\\."))


counts <- read.table(input_args[2], skip=1, sep="\t")
columns <- as.vector(as.matrix(read.table(input_args[2], nrows=1, sep="\t")))
rownames(counts) <- counts[,1]
counts <- counts[,2:(ncol(counts)-1)] 
columns <- columns[1:input_args[5]]

total <- counts[,c(TRUE,FALSE)]
meth <- counts[,c(FALSE,TRUE)]
colnames(total) <- columns
colnames(meth) <- columns


library(bsseq)

BS_unsmoothed <- BSseq(chr = unlist(strsplit(rownames(counts),split = ":"))[c(TRUE,FALSE,FALSE,FALSE)], pos = as.numeric(unlist(strsplit(rownames(counts),split = ":"))[c(FALSE,TRUE,FALSE,FALSE)]), 
               M = meth, Cov = total, sampleNames = columns)
BS_smoothed <- BSmooth(BS_unsmoothed, verbose = TRUE, ns=1) #Default would be ns=70

G1=columns[c(TRUE,FALSE)]
G2=columns[c(FALSE,TRUE)]


BS_fisher<-fisherTests(BS_smoothed, group1=G2, group2=G1, lookup = NULL, returnLookup = FALSE, mc.cores = 1, verbose = TRUE)

print("Saving p-values into files.")
write.table(BS_fisher, file=paste(input_args[4],"_Fisherpvalue.txt",sep=""), row.names=FALSE, col.names=FALSE)




