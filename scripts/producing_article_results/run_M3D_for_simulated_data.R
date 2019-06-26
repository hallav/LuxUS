#This script is for running the M3D tool (Mayo et al., 2014) for multiple simulated data sets 
#(simulated from LuxUS model) at the same time. Multiple data sets are processed 
#at the same time as the M3D tool requires there to be more windows than one.

#Setting up the environment 

.libPaths("/scratch/work/hallav1/R")
Sys.setenv(R_LIBS="/scratch/work/hallav1/R")


#Arguments to the script

input_args <- commandArgs(TRUE)

#The input arguments should be 
#1: The input data file name identifier (with path!)
#2: The number of simulated data files to be used as input.
#3: Output filename identifier (with path!)
#4: Number of cytosines
#5: Number of replicates (total number is double of this)

#Test that the number of arguments is correct
if (length(input_args)!=5){
  stop("The number of input arguments is incorrect.",call.=FALSE)
}

N_data_files=as.numeric(input_args[2])
N_cytosines=as.numeric(input_args[4])
N_replicates=as.numeric(input_args[5])


library(BiSeq)
library(M3D)

#Go through the given input files and save the data into a single BSRaw object


#Objects for storing the window ranges
range_starts <- rep(0,N_data_files*2)
range_ends <- rep(0,N_data_files*2)
range_chromosomes <- rep(0,N_data_files*2) 
  
#Objects for storing the windows into BSRaw objects
window_starts <- rep(0,N_data_files*N_cytosines*2)
window_ends <- rep(0,N_data_files*N_cytosines*2)
window_chromosomes <- rep(0,N_data_files*N_cytosines*2)

totalReads_all <- matrix(NA,nrow=N_data_files*N_cytosines*2,ncol=N_replicates*2)
methReads_all <-matrix(NA,nrow=N_data_files*N_cytosines*2,ncol=N_replicates*2)

ind_i=1

for (d in 0:1){

  for (i in 1:N_data_files){

    print(paste("Loading data set",i,"diff",d,sep=" "))
  
    luxus_simulated <- read.table(paste(input_args[1],i-1,"_diff",d,".txt", sep = ""),skip=1,header = FALSE, sep = "",row.names=1)
    sample_names <- read.table(paste(input_args[1],i-1,"_diff",d,".txt", sep = ""),nrows=1,header = FALSE, sep = "",stringsAsFactors =FALSE)
  
    coordinate_info <- unlist(strsplit(rownames(luxus_simulated),":"))
    coordinates <- as.numeric(coordinate_info[c(FALSE,TRUE,FALSE,FALSE)])
  
    totalReads_window <- as.matrix(luxus_simulated[,c(TRUE,FALSE)])
    methReads_window <- as.matrix(luxus_simulated[,c(FALSE,TRUE)])
  
    totalReads_sorted <- totalReads_window[,c(seq(1,ncol(totalReads_window),by=2),seq(2,ncol(totalReads_window),by=2))]
    methReads_sorted <- methReads_window[,c(seq(1,ncol(methReads_window),by=2),seq(2,ncol(methReads_window),by=2))]
  
    totalReads_all[((ind_i-1)*N_cytosines+1):(ind_i*N_cytosines),] <- totalReads_sorted
    methReads_all[((ind_i-1)*N_cytosines+1):(ind_i*N_cytosines),] <- methReads_sorted
  
    range_starts[ind_i] <- coordinates[1]
    range_ends[ind_i] <- coordinates[N_cytosines]
    range_chromosomes[ind_i] <- paste("chr",ind_i,sep="")
  
    window_starts[((ind_i-1)*N_cytosines+1):(ind_i*N_cytosines)] <-coordinates
    window_ends[((ind_i-1)*N_cytosines+1):(ind_i*N_cytosines)] <-coordinates
    window_chromosomes[((ind_i-1)*N_cytosines+1):(ind_i*N_cytosines)] <- rep(paste("chr",ind_i,sep=""),N_cytosines)

    ind_i=ind_i+1
  }
}

#Construct the BSRaw object

print("Constructing BSRaw object.")

metadata <- list(Sequencer = "Instrument", Year = "2019")

colData_all <- DataFrame(group = c(rep(0,N_replicates),rep(1,N_replicates)),row.names = sample_names[c(seq(1,ncol(totalReads_window),by=2),seq(2,ncol(totalReads_window),by=2))])
colnames(totalReads_all)<-sample_names[c(seq(1,ncol(totalReads_window),by=2),seq(2,ncol(totalReads_window),by=2))]
colnames(methReads_all)<-sample_names[c(seq(1,ncol(methReads_window),by=2),seq(2,ncol(methReads_window),by=2))]
rowRanges_all <- GRanges(seqnames = window_chromosomes,ranges = IRanges(start = window_starts, end = window_ends))

luxus_simulated_BSraw<-BSraw(metadata = metadata,
                             rowRanges = rowRanges_all,
                             colData = colData_all,
                             totalReads = totalReads_all,
                             methReads = methReads_all)

#Run M3D

print("Running M3D.")

test_regions<-GRanges(seqnames = range_chromosomes,ranges = IRanges(start = range_starts, end = range_ends))

overlaps <- findOverlaps(test_regions,rowRanges(luxus_simulated_BSraw))
MMDlist_simulated <- M3D_Wrapper(luxus_simulated_BSraw, overlaps)

head(MMDlist_simulated$Full)

M3Dstat <- MMDlist_simulated$Full-MMDlist_simulated$Coverage

group1 <- unique(colData(luxus_simulated_BSraw)$group)[1]
group2 <-unique(colData(luxus_simulated_BSraw)$group)[2]
pvals_luxus_simulated <- pvals(luxus_simulated_BSraw, test_regions, M3Dstat,group1, group2, smaller=FALSE, comparison='allReps', method='empirical', closePara=0.005)

print("The FDR corrected p-values are:")
head(pvals_luxus_simulated$FDRmean)

print("Saving p-values into files.")
write.table(pvals_luxus_simulated$FDRmean, file=paste(input_args[3],"_diff","both","_FDRmean.txt",sep=""), row.names=FALSE, col.names=FALSE)
write.table(pvals_luxus_simulated$Pmean, file=paste(input_args[3],"_diff","both","_Pmean.txt",sep=""), row.names=FALSE, col.names=FALSE)



