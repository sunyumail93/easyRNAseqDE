#' Merge salmon results to a table
#' @param SampleInfoFile File name of the Sample information text file.
#' @param uniqueMatchingFile Transcript to Gene unique matching file with two columns.
#' @param OutputPrefix Prefix of output files, default is "default".
#' @param write_to_file True of False, write to a file or not. If True, then use OutputPrefix as file prefix.
#' @return Merged salmon Counts and TPM table.
#' @export
#
# This function merges salmon outputs from PipeRNAseq.sh and summarizes the results
#    into a Gene and TPM table
# The output table can be used for further analysis, but not recommended for DE

MergeSalmon <- function(SampleInfoFile, uniqueMatchingFile, OutputPrefix="default", write_to_file=F) {

  print("MergeSalmon Running")

  #Get names
  OutputFileTXT <- paste(OutputPrefix,".salmon.GeneCountsTPM.txt",sep = "")

  #Preparing analysis
  Info <- read.table(SampleInfoFile,col.names = c("Data","DataShortName","condition"),header = F)
  dddn <- dim(Info)[1]
  tx2gene <- read.table(uniqueMatchingFile)
  Files <- as.character(Info$Data)
  names(Files) <- as.character(Info$DataShortName)

  txi <- tximport(Files, type="salmon", tx2gene=tx2gene)
  print("Data resolved:")
  head(txi$counts)

  #Merge tables
  GeneNames <- row.names(txi$abundance)
  Merged <- data.frame(txi$counts,GeneNames,txi$abundance)
  colnames(Merged) <- c(colnames(txi$counts),"GeneNames",colnames(txi$counts))
  if (write_to_file == T){
  write.table(Merged, file=OutputFileTXT,quote = F,col.names = NA, sep="\t")
  }

  print("MergeSalmon Done.")
  return(Merged)
}
