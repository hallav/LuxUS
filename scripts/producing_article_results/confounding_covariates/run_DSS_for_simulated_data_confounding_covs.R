#This script is for running the DSS tool for data sets 
#(simulated from LuxUS model) with confounding covariates.

#Setting up the environment 

.libPaths("/scratch/work/hallav1/R_3_6_1")
Sys.setenv(R_LIBS="/scratch/work/hallav1/R_3_6_1")

#Arguments to the script

input_args <- commandArgs(TRUE)
#The input arguments should be
#1: Folder to be set as environment
#2: Input file name (proportion table)
#3: Input file name (design matrix)
#4: Output file name (+ path if different than environment)
#5: Number of columns in proportion table = total number of replicates


#Test that the number of input arguments is correct
if (length(input_args)!=5) {
  stop("The number of input arguments is incorrect.", call.=FALSE)
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

#library(BiSeq)
#library(M3D)

library(DSS)
require(bsseq)

datalist=list()
for(i in c(1:input_args[5])){
  sample_df=data.frame(chr=rep(1,nrow(counts)), #unlist(strsplit(rownames(counts),split = ":"))[c(TRUE,FALSE,FALSE,FALSE)],
                       pos=as.numeric(unlist(strsplit(rownames(counts),split = ":"))[c(FALSE,TRUE,FALSE,FALSE)]),
                       N=as.matrix(counts[,i*2-1]),X=as.matrix(counts[,i*2]))
  datalist[[i]]=sample_df
}


BSobj=makeBSseqData(datalist,columns)

design_table<- read.table(input_args[3], skip=1, sep="\t")
rownames(design_table) <- design_table[,1]
design_table <- design_table[2:(ncol(design_table))]
IsCase <- as.vector(design_table[,2])
cov1 <- as.vector(design_table[,3])
cov2 <- as.vector(design_table[,4])
DM <- data.frame(IsCase,cov1,cov2) 
#How the intercept term should be handled? -> formula in DMLfit.multiFactor includes intercept
#Source: https://www.bioconductor.org/packages/release/bioc/vignettes/DSS/inst/doc/DSS.html, section 3.4.2

DMLfit = DMLfit.multiFactor(BSobj, design=DM, formula=~IsCase+cov1+cov2)
DMLfit_500bp = DMLfit.multiFactor(BSobj, design=DM, formula=~IsCase+cov1+cov2,smoothing=TRUE,smoothing.span=500)
DMLfit_1000bp = DMLfit.multiFactor(BSobj, design=DM, formula=~IsCase+cov1+cov2,smoothing=TRUE,smoothing.span=1000)
#Index for the covariate to be tested is 2 as intercept is included in the model 
DMLtestresult = DMLtest.multiFactor(DMLfit, coef=2)
DMLtestresult_smoothed_1000bp=DMLtest.multiFactor(DMLfit_1000bp, coef=2) 
DMLtestresult_smoothed_500bp=DMLtest.multiFactor(DMLfit_500bp, coef=2)

#Check that the columns that are being stored are the correct ones!
print("Saving p-values into files.")
write.table(DMLtestresult[,c(1,2,4)], file=paste(input_args[4],"_unsmoothed.txt",sep=""), row.names=FALSE, col.names=FALSE)
write.table(DMLtestresult_smoothed_1000bp[,c(1,2,4)], file=paste(input_args[4],"_smoothed1000bp.txt",sep=""), row.names=FALSE, col.names=FALSE)
write.table(DMLtestresult_smoothed_500bp[,c(1,2,4)], file=paste(input_args[4],"_smoothed500bp.txt",sep=""), row.names=FALSE, col.names=FALSE)


