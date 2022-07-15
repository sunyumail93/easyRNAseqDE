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

  print("Start running easyDE_FromSalmon")

  # Preparing analysis
  print("Reading the SampleInfo file...")
  Info <- read.table(SampleInfo,col.names = c("Data","DataShortName","condition"),header = F)
  dddn <- dim(Info)[1]
  tx2gene <- read.table(uniqueMatchingFile)
  Files <- as.character(Info$Data)
  names(Files) <- as.character(Info$DataShortName)

  # txi$counts are the raw counts, and txi$abundance are the TPM values
  txi <- tximport::tximport(Files, type="salmon", tx2gene=tx2gene)
  print("Data resolved.")
  samples <- Info$DataShortName
  condition <- Info$condition
  coldata <- data.frame(samples, condition)

  # Fitting model
  print("Merging samples...")
  ddsTxi <- DESeq2::DESeqDataSetFromTximport(txi,
                                     colData = coldata,
                                     design = ~ condition)
  dds <- DESeq2::DESeq(ddsTxi)
  rld <- DESeq2::rlog(dds, blind = FALSE)   #This step takes time!

  if (createQuickomicsFiles == T){
    # Exp_data file
    write.csv(SummarizedExperiment::assay(rld), file = paste0(QuickomicsPrefix, "_Exp_rlogData.csv"))
    write.csv(log2(txi$abundance+1), file = paste0(QuickomicsPrefix, "_Exp_Log2TPMData.csv"))
    write.csv(txi$abundance, file = paste0(QuickomicsPrefix, "_Exp_RawTPMData.csv"))

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
  print("Reading the Comparison file...")
  Comparison <- read.table(ComparisonFile,header=F,col.names = c("Condition_1","Condition_2","OutputFileName"))
  print("Going to perform DE analysis on the following comparisons:")
  print(Comparison)

  print("Starting DE analysis")
  print("*****************************************************************")

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

    txi <- tximport::tximport(Files, type="salmon", tx2gene=tx2gene)
    print("Data resolved.")

    #DESeq2 analysis
    samples <- Info_Sub$DataShortName
    condition <- Info_Sub$condition
    coldata <- data.frame(samples, condition)
    print("Running DE analysis...")
    ddsTxi <- DESeq2::DESeqDataSetFromTximport(txi,
                                       colData = coldata,
                                       design = ~ condition)
    dds <- DESeq2::DESeq(ddsTxi)
    res <- DESeq2::results(dds, contrast = c("condition", Condition_2, Condition_1))
    res_ordered <- res[order(res$padj),]
    res_ordered

    CountGenes(res_ordered)
    #Log transformation
    rld <- DESeq2::rlog(dds, blind = FALSE)   #This step takes time!
    #get counts
    counts_normalized <- DESeq2::counts(dds,normalized=TRUE)

    #Combining all results
    padjCutoff <- 0.05
    GeneName <- row.names(txi$counts)
    CombinedAllInfo <- cbind(as.data.frame(res),GeneName,txi$counts,GeneName,txi$abundance,GeneName,as.data.frame(SummarizedExperiment::assay(rld)),GeneName,counts_normalized)
    CombinedAllInfo_sorted <- CombinedAllInfo[order(CombinedAllInfo$padj),]
    OutputFileNameAllCombinedSortedTXT <- paste(OutputPrefix,".salmon.All.DEseq2GeneCountsTPMrlogNormcounts.txt",sep = "")  #This is the final everything combined output
    write.table(CombinedAllInfo_sorted,file = OutputFileNameAllCombinedSortedTXT,quote = F, sep = "\t", col.names = NA)

    print("Start plotting figures...This may take some time!")
    print("Running log transformation")

    #Draw Distance Plot
    colors <- colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
    sampleDists <- dist(t(SummarizedExperiment::assay(rld)))
    sampleDistMatrix <- as.matrix( sampleDists )
    PlotNamesDistance <- paste(OutputPrefix,".salmon.Plot.DistancePlot.All.pdf",sep = "")
    print("Starting figure generation...")
    print("Figure 1: Distance plot")
    pdf(PlotNamesDistance,width = 11,height = 10,onefile=FALSE)
    pheatmap::pheatmap(sampleDistMatrix,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists,
             col = colors)
    dev.off()

    #Draw PCA Plots
    print("Figure 2: PCA")
    plot_pca <- DESeq2::plotPCA(rld)
    ggplot2::ggsave(paste(OutputPrefix,".salmon.Plot.PCA.1.pdf",sep = ""), plot_pca)
    PCA_Data <- DESeq2::plotPCA(rld,intgroup = c("condition", "samples"),returnData = TRUE)
    plot2 <- ggplot2::ggplot(PCA_Data, ggplot2::aes(x = `PC1`, y = `PC2`, color = samples, shape = condition)) +
      ggplot2::geom_point(size = 3) + ggplot2::coord_fixed()
    ggplot2::ggsave(paste(OutputPrefix,".salmon.Plot.PCA.2.pdf",sep = ""), plot2)

    #Draw MA Plot
    print("Figure 3: MA plot")
    PlotNamesMA <- paste(OutputPrefix,".salmon.Plot.MA.pdf",sep = "")
    pdf(PlotNamesMA,width = 12,height = 10)
    res_MA <- tryCatch(expr = {res_MA <- DESeq2::results(dds, addMLE=TRUE)},error=function(e) DESeq2::results(dds))
    MaxValue <- max(res_MA$lfcMLE,na.rm = T) + 0.5
    if (is.finite(MaxValue) == F){
      MaxValue <- max(res_MA$lfcSE[is.finite(res_MA$lfcSE)],na.rm = T) + 0.5
      DESeq2::plotMA(res_MA,main="MA-Plot of shrunken log2 fold changes",ylim=c(-MaxValue,MaxValue))
    }else{
      DESeq2::plotMA(res_MA,MLE=TRUE,main="MA-Plot of shrunken log2 fold changes",ylim=c(-MaxValue,MaxValue))
    }
    TopUp <- row.names(subset(res_MA,log2FoldChange > 0 & padj < 0.05))
    TopDown <- row.names(subset(res_MA,log2FoldChange < 0 & padj < 0.05))
    if (length(TopUp) > 0){
      if (length(TopUp) > 10){
        Top10Up <- TopUp[1:10]
      }else{
        Top10Up <- TopUp
      }
      for (j in 1:dim(res_MA[Top10Up, ])[1]){
        points(res_MA[Top10Up, ][j,1], res_MA[Top10Up, ][j,2], col="red", cex=1, lwd=1)
        text(res_MA[Top10Up, ][j,1], res_MA[Top10Up, ][j,2], Top10Up[j], pos=2, col="red",cex=0.5)
      }
    }
    if (length(TopDown) > 0){
      if (length(TopDown) > 10){
        Top10Down <- TopDown[1:10]
      }else{
        Top10Down <- TopDown
      }
      for (k in 1:dim(res_MA[Top10Down, ])[1]){
        points(res_MA[Top10Down, ][k,1], res_MA[Top10Down, ][k,2], col="red", cex=1, lwd=1)
        text(res_MA[Top10Down, ][k,1], res_MA[Top10Down, ][k,2], Top10Down[k], pos=2, col="dodgerblue",cex=0.5)
      }
    }
    dev.off()

    #Draw volcano Plot
    print("Figure 4: Volcano plot")
    padjCutoff=0.05
    res_ordered_vp <- as.data.frame(res_ordered)
    res_ordered_vp_up <- subset(res_ordered_vp, res_ordered_vp$log2FoldChange > 0 & res_ordered_vp$padj < padjCutoff)
    res_ordered_vp_down <- subset(res_ordered_vp, res_ordered_vp$log2FoldChange < 0 & res_ordered_vp$padj < padjCutoff)
    res_ordered_vp_nochange <- subset(res_ordered_vp, res_ordered_vp$padj >= padjCutoff)
    tryCatch(res_ordered_vp_up$Change <- "Up", error=function(e) {})
    tryCatch(res_ordered_vp_down$Change <- "Down", error=function(e) {})
    tryCatch(res_ordered_vp_nochange$Change <- "No change", error=function(e) {})
    res_ordered_vp_final <- rbind(res_ordered_vp_up, res_ordered_vp_down, res_ordered_vp_nochange)   #The NA values in padj have been removed
    ToPlotList <- c(row.names(res_ordered_vp_up)[1:10], row.names(res_ordered_vp_down)[1:10])
    ToPlotListData <- res_ordered_vp_final[ToPlotList,]
    plot3 <- ggplot2::ggplot(res_ordered_vp_final, ggplot2::aes(x=log2FoldChange, y=-log10(padj), color=Change, size=-log10(padj))) +
      ggplot2::geom_point(alpha=0.4) +
      ggplot2::scale_color_manual(labels = c("Down", "No change", "Up"), values=c("blue", "black", "red")) +
      ggplot2::xlab(expression('Log'['2']*' Fold Change')) + ggplot2::ylab(expression('-Log'['10']*' (adjusted '* italic('p-value') * ')')) +
      ggplot2::theme_classic() +
      ggplot2::geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "grey") +
      ggplot2::geom_vline(xintercept=c(-2,2), linetype="dashed", color = "grey") +
      ggrepel::geom_text_repel(
        data=ToPlotListData,
        ggplot2::aes(label = rownames(ToPlotListData)),
        size          = 4,
        box.padding   = 1.5,
        point.padding = 0.5,
        segment.color = "grey50")
    ggplot2::ggsave(paste(OutputPrefix,".salmon.Plot.VolcanoPlot.pdf",sep = ""), width = 12,height = 10, plot3)

    #Plot a top 50 genes' heatmap
    print("Figure 5: Heatmap")
    df <- data.frame(SummarizedExperiment::colData(dds))
    df <- df[c("condition")]
    heatColors <- colorRampPalette(c("blue", "white", "red"))(n = 500)
    pdf(paste(OutputPrefix,".salmon.Plot.HeatmapTop50.pdf",sep = ""), width=30, height=30, onefile=FALSE)
    pheatmap::pheatmap(SummarizedExperiment::assay(rld)[rownames(res_ordered)[1:50], ], scale="row", show_rownames=TRUE, annotation_col=df, col=heatColors, cellwidth = 50,cellheight = 35,legend = T,
             fontsize_row=25,fontsize_col = 25,fontsize=25)
    dev.off()

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
      }else{
        Comparison_old <- Comparison_current
        Comparison_current <- as.data.frame(res)
        Comparison_current$UniqueID <- row.names(Comparison_current)
        Comparison_current$test <- OutputPrefix
        Comparison_current$Adj.P.Value <- Comparison_current$padj
        Comparison_current$P.Value <- Comparison_current$pvalue
        Comparison_current$logFC <- Comparison_current$log2FoldChange
        Comparison_current <- Comparison_current[,c("UniqueID", "test", "Adj.P.Value", "P.Value", "logFC")]
        Comparison_current <- rbind(Comparison_old, Comparison_current)
      }
    }

    print(paste0("Done running comparison: ",as.character(Compare_Sub$Condition_2)," v.s. ",as.character(Compare_Sub$Condition_1)))
    print("*****************************************************************")
    cat("\n")
  }
  if (createQuickomicsFiles == T){
    write.csv(Comparison_current, file = paste0(QuickomicsPrefix, "_Comparison_data.csv"), row.names = F)
  }
  print("Pipeline finished.")
}
