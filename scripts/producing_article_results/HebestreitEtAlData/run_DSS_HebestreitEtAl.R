#Setting up the environment 

.libPaths("/scratch/work/hallav1/R")
Sys.setenv(R_LIBS="/scratch/work/hallav1/R")

#Arguments to the script

input_args <- commandArgs(TRUE)
#The input arguments should be
#1: Folder to be set as environment
#2: Input file name (proportion table)
#3: Input file name (design matrix) NOT NEEDED
#4: Output file name (+ path if different than environment)
#5: Number of columns in proportion table = total number of replicates


#Test that the number of input arguments is correct
if (length(input_args)!=5) {
  stop("The number of input arguments is incorrect.", call.=FALSE)
}

setwd(input_args[1])

file_id=unlist(strsplit(input_args[2],"\\."))


counts <- read.table(input_args[2], skip=1, sep="")
columns <- as.vector(as.matrix(read.table(input_args[2], nrows=1, sep="\t")))
rownames(counts) <- counts[,1]
counts <- counts[,2:(ncol(counts))] 
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
  sample_df=data.frame(chr=unlist(strsplit(rownames(counts),split = ":"))[c(TRUE,FALSE,FALSE,FALSE)],
                       pos=as.numeric(unlist(strsplit(rownames(counts),split = ":"))[c(FALSE,TRUE,FALSE,FALSE)]),
                       N=as.matrix(counts[,i*2-1]),X=as.matrix(counts[,i*2]))
  datalist[[i]]=sample_df
}


BSobj=makeBSseqData(datalist,columns)
G1=columns[c(rep(TRUE,6),rep(FALSE,6))]
G2=columns[c(rep(FALSE,6),rep(TRUE,6))]

DMLtestresult=DMLtest(BSobj,group1=G1,group2=G2)
DMLtestresult_smoothed_1000bp=DMLtest(BSobj,group1=G1,group2=G2,smoothing=TRUE,smoothing.span=1000)
DMLtestresult_smoothed_2000bp=DMLtest(BSobj,group1=G1,group2=G2,smoothing=TRUE,smoothing.span=2000)
DMLtestresult_smoothed_500bp=DMLtest(BSobj,group1=G1,group2=G2,smoothing=TRUE,smoothing.span=500)

print("Saving p-values into files.")
write.table(DMLtestresult[,c(1,2,10)], file=paste(input_args[4],"_unsmoothed.txt",sep=""), row.names=FALSE, col.names=FALSE)
write.table(DMLtestresult_smoothed_1000bp[,c(1,2,10)], file=paste(input_args[4],"_smoothed1000bp.txt",sep=""), row.names=FALSE, col.names=FALSE)
write.table(DMLtestresult_smoothed_2000bp[,c(1,2,10)], file=paste(input_args[4],"_smoothed2000bp.txt",sep=""), row.names=FALSE, col.names=FALSE)
write.table(DMLtestresult_smoothed_500bp[,c(1,2,10)], file=paste(input_args[4],"_smoothed500bp.txt",sep=""), row.names=FALSE, col.names=FALSE)


