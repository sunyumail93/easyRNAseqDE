#' Run DE analysis using salmon results
#' @param SampleInfo File name of the Sample information text file.
#' @param uniqueMatchingFile Transcript to Gene unique matching file with two columns
#' @param ComparisonFile DE comparison file with three columns: Condition_1 Condition_2 OutputPrefix
#' @param createQuickomicsFiles True or False, output Quickomics files
#' @param QuickomicsPrefix Prefix for Quickomics files
#' @export
#
# This function takes salmon sf output files and two additional files (comparison, and condition annotation)
#    and run the full DESeq2 analysis
# The output results include summary txt files, multiple plots, and Quickomics input files

easyDE_FromSalmon <- function(SampleInfo, uniqueMatchingFile, ComparisonFile, createQuickomicsFiles=F, QuickomicsPrefix=NULL) {

  print("Running easyDE_FromSalmon")

  #Preparing analysis
  Info <- read.table(SampleInfo,col.names = c("Data","DataShortName","condition"),header = F)
  dddn <- dim(Info)[1]
  tx2gene <- read.table(uniqueMatchingFile)
  Files <- as.character(Info$Data)
  names(Files) <- as.character(Info$DataShortName)

  # txi$counts are the raw counts, and txi$abundance are the TPM values
  txi <- tximport(Files, type="salmon", tx2gene=tx2gene)
  print("Data resolved:")
  samples <- Info$DataShortName
  condition <- Info$condition
  coldata <- data.frame(samples, condition)

  if (createQuickomicsFiles == T){
    # Exp_data file
    write.csv(log2(txi$abundance+1), file = paste0(QuickomicsPrefix, "_Exp_Log2Data.csv"))
    write.csv(txi$abundance, file = paste0(QuickomicsPrefix, "_Exp_RawData.csv"))

    # sample meta
    colnames(coldata) <- c("sampleid","group")
    write.csv(coldata, file = paste0(QuickomicsPrefix, "_Sample_metadata.csv"), row.names = F)

    # Protein name optional file
    FinalCounts_temp <- as.data.frame(txi$abundance)
    FinalCounts_temp$id <- 1:dim(FinalCounts_temp)[1]
    FinalCounts_temp$UniqueID <- row.names(FinalCounts_temp)
    FinalCounts_temp$Gene.Name <- row.names(FinalCounts_temp)
    FinalCounts_temp$Protein.ID <- ""
    FinalCounts_temp$GeneType <- "protein_coding"
    Quickomics_optional <- FinalCounts_temp[,c("id", "UniqueID", "Gene.Name", "Protein.ID", "GeneType")]
    write.csv(Quickomics_optional, file = paste0(QuickomicsPrefix, "_ProteinGeneName_optional.csv"), row.names = F)
  }

  #Get comparison
  Comparison <- read.table(ComparisonFile,header=F,col.names = c("Condition_1","Condition_2","OutputFileName"))
  Comparison

  print("*****************************************************************")
  print("Starting DE analysis")

  #Paired differential expression analysis
  for (i in 1:dim(Comparison)[1]){
    Compare_Sub <- Comparison[i,]
    print(paste0("Reading the paired comparison ",i,": ",as.character(Compare_Sub$Condition_2)," v.s. ",as.character(Compare_Sub$Condition_1)))
    Condition_1 <- as.character(Comparison[i,1])
    Condition_2 <- as.character(Comparison[i,2])
    OutputPrefix <- as.character(Comparison[i,3])
    LabelCondition_1 <- subset(Info,Info$condition==Condition_1)
    LabelCondition_2 <- subset(Info,Info$condition==Condition_2)
    OriDataLabel <- c(as.character(LabelCondition_1$Data),as.character(LabelCondition_2$Data))
    samples <- c(as.character(LabelCondition_1$DataShortName),as.character(LabelCondition_2$DataShortName))
    condition <- c(as.character(LabelCondition_1$condition),as.character(LabelCondition_2$condition))

    #Create Info_Sub
    Info_Sub <- data.frame(OriDataLabel, samples, condition)
    colnames(Info_Sub) <- c("Data", "DataShortName", "condition")
    Files <- as.character(Info_Sub$Data)
    names(Files) <- as.character(Info_Sub$DataShortName)

    txi <- tximport(Files, type="salmon", tx2gene=tx2gene)
    print("Data resolved:")

    #DESeq2 analysis
    samples <- Info_Sub$DataShortName
    condition <- Info_Sub$condition
    coldata <- data.frame(samples, condition)
    ddsTxi <- DESeqDataSetFromTximport(txi,
                                       colData = coldata,
                                       design = ~ condition)
    dds <- DESeq(ddsTxi)
    res <- results(dds, contrast = c("condition", Condition_2, Condition_1))
    res_ordered <- res[order(res$padj),]
    res_ordered

    CountGenes(res_ordered)

    #Combining all results
    padjCutoff <- 0.05
    GeneName <- row.names(txi$counts)
    CombinedAllInfo <- cbind(as.data.frame(res),GeneName,txi$counts,GeneName,txi$abundance)
    CombinedAllInfo_sorted <- CombinedAllInfo[order(CombinedAllInfo$padj),]
    OutputFileNameAllCombinedSortedTXT <- paste(OutputPrefix,".salmon.All.DEseq2GeneCountsTPM.txt",sep = "")  #This is the final everything combined output
    write.table(CombinedAllInfo_sorted,file = OutputFileNameAllCombinedSortedTXT,quote = F, sep = "\t", col.names = NA)

    print("Start plotting figures...This may take some time!")
    print("Running log transformation")
    #Log transformation
    rld <- rlog(dds, blind = FALSE)   #This step takes time!

    #Draw Distance Plot
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    sampleDists <- dist(t(assay(rld)))
    sampleDistMatrix <- as.matrix( sampleDists )
    PlotNamesDistance <- paste(OutputPrefix,".salmon.Plot.DistancePlot.All.pdf",sep = "")
    pdf(PlotNamesDistance,width = 11,height = 10,onefile=FALSE)
    pheatmap(sampleDistMatrix,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists,
             col = colors)
    dev.off()

    #Draw MA Plot
    print("Plotting MA")
    PlotNamesMA <- paste(OutputPrefix,".salmon.Plot.MA.pdf",sep = "")
    pdf(PlotNamesMA,width = 12,height = 10)
    res_MA <- tryCatch(expr = {res_MA <- results(dds, addMLE=TRUE)},error=function(e) results(dds))
    MaxValue <- max(res_MA$lfcMLE,na.rm = T) + 0.5
    if (is.finite(MaxValue) == F){
      MaxValue <- max(res_MA$lfcSE[is.finite(res_MA$lfcSE)],na.rm = T) + 0.5
      plotMA(res_MA,main="MA-Plot of shrunken log2 fold changes",ylim=c(-MaxValue,MaxValue))
    }else{
      plotMA(res_MA,MLE=TRUE,main="MA-Plot of shrunken log2 fold changes",ylim=c(-MaxValue,MaxValue))
    }
    dev.off()

    #Draw PCA Plots
    print("Plotting PCA")
    plot_pca <- plotPCA(rld)
    ggsave(paste(OutputPrefix,".salmon.Plot.PCA.1.pdf",sep = ""), plot_pca)
    PCA_Data <- plotPCA(rld,intgroup = c("condition", "samples"),returnData = TRUE)
    plot2 <- ggplot(PCA_Data, aes(x = `PC1`, y = `PC2`, color = samples, shape = condition)) +
      geom_point(size = 3) + coord_fixed()
    ggsave(paste(OutputPrefix,".salmon.Plot.PCA.2.pdf",sep = ""), plot2)

    #Draw volcano Plot
    print("Plotting Volcano, using padj<0.05 as significant")
    PlotNamesVolcano <- paste(OutputPrefix,".salmon.Plot.Volcano.pdf",sep = "")
    res_ordered_vp <- as.data.frame(res_ordered)
    res_ordered_vp$threshold <- as.factor( res_ordered_vp$padj < padjCutoff) #You can adjust threshold here
    res_ordered_vp_noNA <- subset(res_ordered_vp,threshold != 'NA')
    volcano <- ggplot(res_ordered_vp_noNA, aes(x=log2FoldChange, y=-log10(pvalue), color=threshold)) +
      geom_point(alpha=0.4, size=0.8) +
      xlab("log2 fold change") + ylab("-log10 p-value")
    ggsave(paste(OutputPrefix,".salmon.Plot.Volcano.pdf",sep = ""), volcano)

    #Create Quickomics files
    if (createQuickomicsFiles == T){
      if (i == 1){
        Comparison_current <- as.data.frame(res)
        Comparison_current$UniqueID <- row.names(Comparison_current)
        Comparison_current$test <- OutputPrefix
        Comparison_current$Adj.P.Value <- Comparison_current$padj
        Comparison_current$P.Value <- Comparison_current$pvalue
        Comparison_current$logFC <- Comparison_current$log2FoldChange
        Comparison_current <- Comparison_current[,c("UniqueID", "test", "Adj.P.Value", "P.Value", "logFC")]
        print(i)
        print(dim(Comparison_current))
      }else{
        print(i)
        Comparison_old <- Comparison_current
        Comparison_current <- as.data.frame(res)
        Comparison_current$UniqueID <- row.names(Comparison_current)
        Comparison_current$test <- OutputPrefix
        Comparison_current$Adj.P.Value <- Comparison_current$padj
        Comparison_current$P.Value <- Comparison_current$pvalue
        Comparison_current$logFC <- Comparison_current$log2FoldChange
        Comparison_current <- Comparison_current[,c("UniqueID", "test", "Adj.P.Value", "P.Value", "logFC")]
        Comparison_current <- rbind(Comparison_old, Comparison_current)
        print(dim(Comparison_current))
      }
    }

    print(paste0("Done running comparisons: ",as.character(Compare_Sub$Condition_2)," v.s. ",as.character(Compare_Sub$Condition_1)))
    print("*****************************************************************")
    cat("\n\n")
  }
  if (createQuickomicsFiles == T){
    write.csv(Comparison_current, file = paste0(QuickomicsPrefix, "_Comparison_data.csv"), row.names = F)
  }
}
