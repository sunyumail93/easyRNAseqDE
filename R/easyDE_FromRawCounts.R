#' Run DE analysis using raw counts
#' @param count_matrix data.frame of Raw count matrix
#' @param LabelFile File name of the Sample labels, with three columns: Data name, Data short name, condition
#' @param ComparisonFile DE comparison file with three columns: Condition_1 Condition_2 OutputPrefix
#' @param createQuickomicsFiles True or False, output Quickomics files
#' @param QuickomicsPrefix Prefix for Quickomics files
#' @export
#
# This function takes raw count table and two additional files (comparison, and condition annotation)
#    and run the full DESeq2 analysis
# The output results include summary txt files, multiple plots, and Quickomics input files

easyDE_FromRawCounts <- function(count_matrix, LabelFile, ComparisonFile, createQuickomicsFiles=F, QuickomicsPrefix=NULL) {

  print("Running easyDE_FromRawCounts")
  #count_matrix <- as.data.frame(read.table(CountTableFile,header=TRUE,check.names = F))
  Label <- as.data.frame(read.table(LabelFile,header=F,col.names = c("Data","DataShortName","Condition")))
  print("Reading Label file...")
  Label
  Comparison <- read.table(ComparisonFile,header=F,col.names = c("Condition_1","Condition_2","OutputFileName"))
  print("Reading Comparison file...")
  print(Comparison)

  #Prepare data table
  FinalCounts <- count_matrix[,as.character(Label$Data)]
  colnames(FinalCounts) <- Label$DataShortName
  coldata <- Label[,c("DataShortName","Condition")]
  colnames(coldata) <- c("samples","condition")
  condition <- as.character(coldata$condition)

  if (createQuickomicsFiles == T){
  # Exp_data file
  write.csv(log2(FinalCounts+1), file = paste0(QuickomicsPrefix, "_Exp_Log2Data.csv"))
  write.csv(FinalCounts, file = paste0(QuickomicsPrefix, "_Exp_RawData.csv"))

  # sample meta
  colnames(coldata) <- c("sampleid","group")
  write.csv(coldata, file = paste0(QuickomicsPrefix, "_Sample_metadata.csv"), row.names = F)

  # Protein name optional file
  FinalCounts_temp <- FinalCounts
  FinalCounts_temp$id <- 1:dim(FinalCounts_temp)[1]
  FinalCounts_temp$UniqueID <- row.names(FinalCounts_temp)
  FinalCounts_temp$Gene.Name <- row.names(FinalCounts_temp)
  FinalCounts_temp$Protein.ID <- ""
  FinalCounts_temp$GeneType <- "protein_coding"
  Quickomics_optional <- FinalCounts_temp[,c("id", "UniqueID", "Gene.Name", "Protein.ID", "GeneType")]
  write.csv(Quickomics_optional, file = paste0(QuickomicsPrefix, "_ProteinGeneName_optional.csv"), row.names = F)
  }

  print("*****************************************************************")
  print("Starting DE analysis")

  #Paired differential expression analysis
  for (i in 1:dim(Comparison)[1]){
    Compare_Sub <- Comparison[i,]
    print(paste0("Reading the paired comparison ",i,": ",as.character(Compare_Sub$Condition_2)," v.s. ",as.character(Compare_Sub$Condition_1)))
    Condition_1 <- as.character(Comparison[i,1])
    Condition_2 <- as.character(Comparison[i,2])
    OutputFileName <- as.character(Comparison[i,3])
    LabelCondition_1 <- subset(Label,Label$Condition==Condition_1)
    LabelCondition_2 <- subset(Label,Label$Condition==Condition_2)
    OriDataLabel <- c(as.character(LabelCondition_1$Data),as.character(LabelCondition_2$Data))
    samples <- c(as.character(LabelCondition_1$DataShortName),as.character(LabelCondition_2$DataShortName))

    FinalCounts <- count_matrix[,OriDataLabel]
    colnames(FinalCounts) <- samples
    condition <- c(as.character(LabelCondition_1$Condition),as.character(LabelCondition_2$Condition))
    coldata <- data.frame(samples, condition)
    coldata

    #DESeq2 analysis
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = FinalCounts, colData = coldata, design = ~ condition)
    dds
    dds <- DESeq2::DESeq(dds)
    sizeFactor <- DESeq2::sizeFactors(dds)
    sizeFactor
    res <- DESeq2::results(dds, contrast = c("condition", Condition_2, Condition_1))
    res_ordered <- res[order(res$padj),]
    res_ordered

    CountGenes(res_ordered)

    #Finalize and output the results
    #The FinalOuput has three parts: Original DESeq2 results, Raw counts, Normalized counts, separated by GeneList, sorted by padj
    DESeq2Result <- as.data.frame(res_ordered)
    GeneOrder <- row.names(DESeq2Result)
    FinalCounts_ordered <- FinalCounts[GeneOrder,]
    FinalCounts_ordered_normalized <- t(t(FinalCounts_ordered)/sizeFactor)
    FinalOutput <- cbind(DESeq2Result,GeneOrder,FinalCounts_ordered,GeneOrder,FinalCounts_ordered_normalized)
    head(FinalOutput,n=3)
    FinalOutput_FileName <- paste(OutputFileName,".DESeq2.txt",sep="")
    write.table(FinalOutput,file = FinalOutput_FileName,quote = F,sep = "\t", col.names=NA)

    #Plot1: Distance plot to show the similarity of the data
    #Log transformation
    print("Running log transformation")
    rld <- DESeq2::rlog(dds, blind = FALSE)   #This step may take time.
    head(SummarizedExperiment::assay(rld), 3)
    colors <- colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
    sampleDists <- dist(t(SummarizedExperiment::assay(rld)))
    sampleDistMatrix <- as.matrix( sampleDists )
    print("Starting figure generation...")
    print("Figure 1: Distance plot")
    pdf(paste(OutputFileName,".DistancePlot.pdf",sep = ""),onefile=FALSE)
    pheatmap::pheatmap(sampleDistMatrix,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists,
             col = colors,main="Distance plot using all genes")
    dev.off()

    #Plot2: PCA
    print("Figure 2: PCA")
    plot_pca <- DESeq2::plotPCA(rld)
    ggplot2::ggsave(paste(OutputFileName,".PCA.1.pdf",sep = ""), plot_pca)
    PCA_Data <- DESeq2::plotPCA(rld,intgroup = c( "condition", "samples"),returnData = TRUE)
    plot2 <- ggplot2::ggplot(PCA_Data, ggplot2::aes(x = `PC1`, y = `PC2`, color = samples, shape = condition)) +
      ggplot2::geom_point(size = 3) + ggplot2::coord_fixed() +
      ggplot2::ggtitle("PCA plot using all genes") + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    ggplot2::ggsave(paste(OutputFileName,".PCA.2.pdf",sep = ""), plot2)

    #Plot3: MA plot, do top 10 up and top 10 down significant genes
    print("Figure 3: MA plot")
    res_MA <- tryCatch(expr = {res_MA <- DESeq2::results(dds, addMLE=TRUE, contrast = c("condition",Condition_2,Condition_1))},error=function(e) DESeq2::results(dds, contrast = c("condition",Condition_2,Condition_1)))
    pdf(paste(OutputFileName,".MAPlot.Labeltop20.pdf",sep = ""),width = 12,height = 10)
    MaxValue <- tryCatch(expr = {MaxValue <- max(res_MA$lfcMLE,na.rm = T) + 0.5},error=function(e) MaxValue <- max(res_MA$lfcSE,na.rm = T) + 0.5)
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
      for (i in 1:dim(res_MA[Top10Up, ])[1]){
        points(res_MA[Top10Up, ][i,1], res_MA[Top10Up, ][i,2], col="red", cex=1, lwd=1)
        text(res_MA[Top10Up, ][i,1], res_MA[Top10Up, ][i,2], Top10Up[i], pos=2, col="red",cex=0.5)
      }
    }
    if (length(TopDown) > 0){
      if (length(TopDown) > 10){
        Top10Down <- TopDown[1:10]
      }else{
        Top10Down <- TopDown
      }
      for (i in 1:dim(res_MA[Top10Down, ])[1]){
        points(res_MA[Top10Down, ][i,1], res_MA[Top10Down, ][i,2], col="red", cex=1, lwd=1)
        text(res_MA[Top10Down, ][i,1], res_MA[Top10Down, ][i,2], Top10Down[i], pos=2, col="red",cex=0.5)
      }
    }
    dev.off()

    #Plot4: Volcano plot, use cutoff: padj<0.05
    print("Figure 4: Volcano plot")
    padjCutoff=0.05
    res_ordered_vp <- as.data.frame(res_ordered)
    res_ordered_vp_up <- subset(res_ordered_vp, res_ordered_vp$log2FoldChange > 0 & res_ordered_vp$padj < padjCutoff)
    res_ordered_vp_down <- subset(res_ordered_vp, res_ordered_vp$log2FoldChange < 0 & res_ordered_vp$padj < padjCutoff)
    res_ordered_vp_nochange <- subset(res_ordered_vp, res_ordered_vp$padj >= padjCutoff)
    res_ordered_vp_up$Change <- "Up"
    res_ordered_vp_down$Change <- "Down"
    res_ordered_vp_nochange$Change <- "No change"
    res_ordered_vp_final <- rbind(res_ordered_vp_up, res_ordered_vp_down, res_ordered_vp_nochange)   #The NA values in padj have been removed
    ToPlotList <- c(row.names(res_ordered_vp_up)[1:10], row.names(res_ordered_vp_down)[1:10])
    ToPlotListData <- res_ordered_vp_final[ToPlotList,]
    # pdf(paste(OutputFileName,".VolcanoPlot.pdf",sep = ""),width = 12,height = 10)
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
    ggplot2::ggsave(paste(OutputFileName,".VolcanoPlot.pdf",sep = ""), width = 12,height = 10, plot3)

    #Draw a top 50 genes heatmap
    print("Figure 5: Heatmap")
    df <- data.frame(SummarizedExperiment::colData(dds))
    df <- df[c("condition","sizeFactor")]
    heatColors <- colorRampPalette(c("blue", "white", "red"))(n = 500)
    pdf(paste(OutputFileName,".HeatmapTop50.pdf",sep = ""), width=30, height=30, onefile=FALSE)
    pheatmap::pheatmap(SummarizedExperiment::assay(rld)[rownames(res_ordered)[1:50], ], scale="row", show_rownames=TRUE, annotation_col=df, col=heatColors, cellwidth = 50,cellheight = 35,legend = T,
             fontsize_row=25,fontsize_col = 25,fontsize=25)
    dev.off()

    #Create Quickomics files
    if (createQuickomicsFiles == T){
      if (i == 1){
      Comparison_current <- DESeq2Result
      Comparison_current$UniqueID <- row.names(Comparison_current)
      Comparison_current$test <- OutputFileName
      Comparison_current$Adj.P.Value <- Comparison_current$padj
      Comparison_current$P.Value <- Comparison_current$pvalue
      Comparison_current$logFC <- Comparison_current$log2FoldChange
      Comparison_current <- Comparison_current[,c("UniqueID", "test", "Adj.P.Value", "P.Value", "logFC")]
      }else{
        Comparison_old <- Comparison_current
        Comparison_current <- DESeq2Result
        Comparison_current$UniqueID <- row.names(Comparison_current)
        Comparison_current$test <- OutputFileName
        Comparison_current$Adj.P.Value <- Comparison_current$padj
        Comparison_current$P.Value <- Comparison_current$pvalue
        Comparison_current$logFC <- Comparison_current$log2FoldChange
        Comparison_current <- Comparison_current[,c("UniqueID", "test", "Adj.P.Value", "P.Value", "logFC")]
        Comparison_current <- rbind(Comparison_old, Comparison_current)
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
