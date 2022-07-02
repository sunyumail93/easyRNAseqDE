#' Merge featureCounts results to a table
#' @param SampleInfoFile File name of the Sample information text file.
#' @param OutputPrefix Prefix of output files, default is "default".
#' @param write_to_file True of False, write to a file or not. If True, then use OutputPrefix as file prefix.
#' @return Merged featureCounts table.
#' @export
#
# This function merges featureCounts outputs from PipeRNAseq.sh and summarizes the results
#    into a text table
# The output table can be used for downstream DE analysis

MergeFeatureCounts <- function(SampleInfoFile, OutputPrefix="default", write_to_file=F) {

  print("MergeFeatureCounts Running")
  #Preparing names
  CombinedNameTXT <- paste(OutputPrefix,".FC.combined.txt",sep = "")
  CombinedNameTXTWithLen <- paste(OutputPrefix,".FC.Length.combined.txt",sep = "")
  CombinedNameRPKMTXT <- paste(OutputPrefix,".FC.RPKM.combined.txt",sep = "")

  #Processing
  Info <- read.table(SampleInfoFile,col.names = c("Data","DataShortName","condition"),header = F)
  Files <- as.character(Info$Data)
  names(Files) <- as.character(Info$DataShortName)

  #Read the first data in
  FirstData <- read.table(Files[1],header=T)
  GeneNum <- dim(FirstData)[1]
  GeneNames <- as.character(FirstData[,1])
  FinalCounts <- data.frame(row.names = GeneNames)

  for (i in 1:length(Files)){
    print(paste("Reading data: ",i,sep = ""))
    CurrentData <- read.table(Files[i],header=T)
    CurrentCounts <- CurrentData[,3]
    FinalCounts <- data.frame(FinalCounts,CurrentCounts)
  }
  colnames(FinalCounts) <- names(Files)
  print("Data resolved:")
  head(FinalCounts)

  #Calculate RPKM
  print("Calculating RPKM...")
  Length <- FirstData$Length
  FinalCountsWithLength <- cbind(Length,FinalCounts)
  TotalReads <- apply(FinalCountsWithLength,2,FUN = sum)
  FinalCountsNormDepth <- t(t(FinalCountsWithLength)/TotalReads*1000000)
  FinalCountsNormDepthNormLength <- FinalCountsNormDepth/Length*1000
  FinalCountsNormDepthNormLength_Output <- as.data.frame(FinalCountsNormDepthNormLength[,2:ncol(FinalCountsNormDepthNormLength)])

  #Output
  if (write_to_file == T){
  write.table(FinalCounts,file = CombinedNameTXT,quote = F,sep = "\t", col.names=NA)
  print("Output Raw Table: Done")
  write.table(FinalCountsWithLength,file = CombinedNameTXTWithLen,quote = F,sep = "\t", col.names=NA)
  write.table(FinalCountsNormDepthNormLength_Output,file = CombinedNameRPKMTXT,quote = F,sep = "\t", col.names=NA)
  print("Output RPKM Table: Done")
  }

  print("MergeFeatureCounts Done.")
  return(FinalCounts)
}
