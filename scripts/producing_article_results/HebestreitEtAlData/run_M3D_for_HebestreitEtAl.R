#This script is for running the M3D tool (Mayo et al., 2014) for multiple simulated DMRs 
#(from Hebestreit et al., 2013) at the same time. Multiple data sets are processed 
#at the same time, as the M3D tool requires there to be more windows than one.

#Setting up the environment 

.libPaths("/scratch/work/hallav1/R")
Sys.setenv(R_LIBS="/scratch/work/hallav1/R")


#Arguments to the script

input_args <- commandArgs(TRUE)

#The input arguments should be 
#1: The input data file name identifier (with path!) for DMR
#2: The input data file name identifier (with path!) for nonDMR
#3: The number of simulated data files to be used as input.
#4: Output filename identifier (with path!)
#5: Maximum number of cytosines per window


#Test that the number of arguments is correct
if (length(input_args)!=5){
  stop("The number of input arguments is incorrect.",call.=FALSE)
}

N_data_files=as.numeric(input_args[3])
N_max_cytosines=as.numeric(input_args[5])
N_replicates=6

library(BiSeq)
library(M3D)

#Go through the given input files and save the data into a single BSRaw object


#Objects for storing the window ranges
range_starts <- rep(0,N_data_files*2)
range_ends <- rep(0,N_data_files*2)
range_chromosomes <- rep(0,N_data_files*2) 
  
#Objects for storing the windows into BSRaw objects
window_starts <- rep(0,N_data_files*N_max_cytosines*2)
window_ends <- rep(0,N_data_files*N_max_cytosines*2)
window_chromosomes <- rep(0,N_data_files*N_max_cytosines*2)

totalReads_all <- matrix(NA,nrow=N_data_files*N_max_cytosines*2,ncol=N_replicates*2)
methReads_all <-matrix(NA,nrow=N_data_files*N_max_cytosines*2,ncol=N_replicates*2)

dmr_or_not <- rep(0,N_data_files*2)

ind_i=1
ind_for_chr=1

N_cytosines_until_now <- 1 #This is an indexing for matrices, so it has to start from 1.

for (d in 0:1){

  if (d==0) {

    inputfileformat=input_args[2]

  } else {

    inputfileformat=input_args[1]
  }

  #Purkkaratkaisu
  for (i in 1:N_data_files){

    print(paste("Loading data set",i,"diff",d,sep=" "))
  
    luxus_simulated <- read.table(paste(inputfileformat,i,".txt", sep = ""),skip=1,header = FALSE, sep = "",row.names=1)
    sample_names <- read.table(paste(inputfileformat,i,".txt", sep = ""),nrows=1,header = FALSE, sep = "",stringsAsFactors =FALSE)
  
    coordinate_info <- unlist(strsplit(rownames(luxus_simulated),":"))
    coordinates <- as.numeric(coordinate_info[c(FALSE,TRUE,FALSE,FALSE)])
  
    totalReads_window <- as.matrix(luxus_simulated[,c(TRUE,FALSE)])
    methReads_window <- as.matrix(luxus_simulated[,c(FALSE,TRUE)])

    N_cytosines_window <- dim(totalReads_window)[1]

    #There should be no windows with only one cytosine, so this restriction is not actually needed
    if (N_cytosines_window>1){ 
  
      #Sorting not needed as the samples are already sorted

      totalReads_all[N_cytosines_until_now:(N_cytosines_until_now+N_cytosines_window-1),] <- totalReads_window
      methReads_all[N_cytosines_until_now:(N_cytosines_until_now+N_cytosines_window-1),] <- methReads_window

      range_starts[ind_i] <- coordinates[1]
      range_ends[ind_i] <- coordinates[N_cytosines_window]
      range_chromosomes[ind_i] <- coordinate_info[1]
  
      window_starts[N_cytosines_until_now:(N_cytosines_until_now+N_cytosines_window-1)] <-coordinates
      window_ends[N_cytosines_until_now:(N_cytosines_until_now+N_cytosines_window-1)] <-coordinates
      window_chromosomes[N_cytosines_until_now:(N_cytosines_until_now+N_cytosines_window-1)] <- coordinate_info[1]

      dmr_or_not[ind_i] <- d

      N_cytosines_until_now <- N_cytosines_until_now+N_cytosines_window
      ind_i=ind_i+1


    } else {

      print("There was only one cytosine in the window.")

    }    


    ind_for_chr=ind_for_chr+1

  }
}

#Construct the BSRaw object

print("Constructing BSRaw object.")

metadata <- list(Sequencer = "Instrument", Year = "2019")

colData_all <- DataFrame(group = c(rep(1,N_replicates),rep(0,N_replicates)),row.names = sample_names) #changed to match Klein & Hebestreit
colnames(totalReads_all) <- sample_names
colnames(methReads_all) <- sample_names
rowRanges_all <- GRanges(seqnames = window_chromosomes[1:N_cytosines_until_now-1],ranges = IRanges(start = window_starts[1:N_cytosines_until_now-1], end = window_ends[1:N_cytosines_until_now-1]))

luxus_simulated_BSraw_notSorted<-BSraw(metadata = metadata,
                             rowRanges = rowRanges_all,
                             colData = colData_all,
                             totalReads = totalReads_all[1:N_cytosines_until_now-1,],
                             methReads = methReads_all[1:N_cytosines_until_now-1,])


luxus_simulated_BSraw <- sort(luxus_simulated_BSraw_notSorted)


#Run M3D

print("Running M3D.")

print(range_chromosomes[1:(ind_i-1)])
print(range_starts[1:(ind_i-1)])
print(range_ends[1:(ind_i-1)])
print(range_ends[1:(ind_i-1)]-range_starts[1:(ind_i-1)])


test_regions<-GRanges(seqnames = range_chromosomes[1:(ind_i-1)],ranges = IRanges(start = range_starts[1:(ind_i-1)], end = range_ends[1:(ind_i-1)]))

overlaps <- findOverlaps(test_regions,rowRanges(luxus_simulated_BSraw))

MMDlist_simulated <- M3D_Wrapper(luxus_simulated_BSraw, overlaps)

warnings()

head(MMDlist_simulated$Full)

M3Dstat <- MMDlist_simulated$Full-MMDlist_simulated$Coverage

group1 <- unique(colData(luxus_simulated_BSraw)$group)[1]
group2 <-unique(colData(luxus_simulated_BSraw)$group)[2]
pvals_luxus_simulated <- pvals(luxus_simulated_BSraw, test_regions, M3Dstat,group1, group2, smaller=FALSE, comparison='allReps', method='empirical', closePara=0.005)

print("The FDR corrected p-values are:")
head(pvals_luxus_simulated$FDRmean)

print("Saving p-values into files.")
write.table(pvals_luxus_simulated$FDRmean, file=paste(input_args[4],"_FDRmean_filteredData.txt",sep=""), row.names=FALSE, col.names=FALSE)
write.table(pvals_luxus_simulated$Pmean, file=paste(input_args[4],"_Pmean_filteredData.txt",sep=""), row.names=FALSE, col.names=FALSE)
write.table(dmr_or_not[1:(ind_i-1)], file=paste(input_args[4],"_dmr_or_not.txt",sep=""), row.names=FALSE, col.names=FALSE)


