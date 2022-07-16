#' Run GO and GSEA analysis using DE results.
#' Please run this function after running easyDE_FromRawCounts or easyDE_FromSalmon.
#' If you would like to prepare your own files, the input format should be DESeq2 output table.
#' The output results include GO plots (BP, CC, and MF), GSEA results and plots (Top 10 GSEA).
#' @param ComparisonFile DE comparison file with three columns: Condition_1 Condition_2 OutputPrefix.
#' @param OrganDatabase Organism GO annotation database, for example, org.Mm.eg.db for mouse.
#' @param GOpajd_cutoff pajd cutoff to define significant genes, default=0.05.
#' @param GOlog2FC_cutoff log2 fold change cutoff to define significant genes, default=0.05
#' @param ShowTermNum Number of GOs to plot, default=20.
#' @param forcerun Force to run the analysis if DE gene number<10, default=False.
#' @importFrom enrichplot barplot
#' @export
#
# Run GO and GSEA analysis using DE results. Please run this function after running easyDE_FromRawCounts or easyDE_FromSalmon.
# The output results include GO plots (BP, CC, and MF), GSEA results and plots (Top 10 GSEA)

easyEnrich_GO_GSEA <- function(ComparisonFile, OrganDatabase, GOpajd_cutoff=0.05, GOlog2FC_cutoff=1.5, ShowTermNum=20, forcerun=F){
  print("Start running easyGO_GSEA")
  Comparison <- read.table(ComparisonFile,header=F,col.names = c("Condition_1","Condition_2","OutputFileName"))
  print("*****************************************************************")

  for (i in 2:dim(Comparison)[1]){
    Compare_Sub <- Comparison[i,]
    print(paste0("Working on comparison ",i,": ",as.character(Compare_Sub$Condition_2)," v.s. ",as.character(Compare_Sub$Condition_1)))
    Condition_1 <- as.character(Comparison[i,1])
    Condition_2 <- as.character(Comparison[i,2])
    OutputFileName <- as.character(Comparison[i,3])

    if (file.exists(paste0(OutputFileName, ".salmon.All.DEseq2GeneCountsTPMrlogNormcounts.txt"))){
      Data <- read.table(paste0(OutputFileName, ".salmon.All.DEseq2GeneCountsTPMrlogNormcounts.txt"), sep = "\t", header = T, row.names = 1, check.names = F)
    }else if (file.exists(paste0(OutputFileName, ".RawCounts.All.DEseq2GeneCountsNormcounts.txt"))){
      Data <- read.table(paste0(OutputFileName, ".RawCounts.All.DEseq2GeneCountsNormcounts.txt"), sep = "\t", header = T, row.names = 1, check.names = F)
    }else{
      print("Error: no DE results can be found. Please run easyDE_FromRawCounts or easyDE_FromSalmon first.")
    }

    print("Filtering genes based on cutoffs: ")
    print(paste0("log2 Fold change cutoff: ", GOlog2FC_cutoff))
    print(paste0("padj cutoff: ", GOpajd_cutoff))
    Data_filtered <- subset(Data, abs(Data$log2FoldChange) > GOlog2FC_cutoff & Data$padj < GOpajd_cutoff)[,1:7]
    print(paste0("There are ",dim(Data_filtered)[1]," genes passed the cutoff."))

    print("Converting gene symbol to Entrez id...")
    GeneList_Matched = clusterProfiler::bitr(as.character(Data_filtered$GeneName), fromType="SYMBOL", toType="ENTREZID", OrgDb=OrganDatabase)
    print("The convertion is done. The first a few rows are:")
    print(head(GeneList_Matched))
    if (dim(GeneList_Matched)[1] < 10){
      print(paste0("Only ",dim(GeneList_Matched)[1]," genes have been converted. The results may not be useful!"))
      if (forcerun == T){
        print("Gene number <10, but force to perform analysis on this comparison")
      }else{
        print("Skip this comparison")
        print("*****************************************************************")
        cat("\n")
        next
      }
    }else{
      print(paste0("Total ",dim(GeneList_Matched)[1]," genes have been converted."))
    }

    print("Running GO analysis on: BP")
    GeneList_Enrich_BP <- clusterProfiler::enrichGO(gene          = GeneList_Matched$ENTREZID,
                                                    OrgDb         = OrganDatabase,
                                                    ont           = "BP",      #Can be BP, CC, MF
                                                    pAdjustMethod = "BH",
                                                    pvalueCutoff  = 1,
                                                    qvalueCutoff  = 1,
                                                    readable      = TRUE)
    print("Running GO analysis on: CC")
    GeneList_Enrich_CC <- clusterProfiler::enrichGO(gene          = GeneList_Matched$ENTREZID,
                                                    OrgDb         = OrganDatabase,
                                                    ont           = "CC",      #Can be BP, CC, MF
                                                    pAdjustMethod = "BH",
                                                    pvalueCutoff  = 1,
                                                    qvalueCutoff  = 1,
                                                    readable      = TRUE)
    print("Running GO analysis on: MF")
    GeneList_Enrich_MF <- clusterProfiler::enrichGO(gene          = GeneList_Matched$ENTREZID,
                                                    OrgDb         = OrganDatabase,
                                                    ont           = "MF",      #Can be BP, CC, MF
                                                    pAdjustMethod = "BH",
                                                    pvalueCutoff  = 1,
                                                    qvalueCutoff  = 1,
                                                    readable      = TRUE)
    geneList <- Data_filtered[GeneList_Matched$SYMBOL,2]
    names(geneList) <- GeneList_Matched$SYMBOL
    print("Plotting figures for: BP")
    tryCatch(BP1 <- barplot(GeneList_Enrich_BP,showCategory = ShowTermNum), error=function(e) print("Not able to plot"))
    tryCatch(BP2 <- clusterProfiler::dotplot(GeneList_Enrich_BP,showCategory = ShowTermNum), error=function(e) print("Not able to plot"))
    tryCatch(BP3 <- clusterProfiler::cnetplot(GeneList_Enrich_BP), error=function(e) print("Not able to plot"))
    tryCatch(BP4 <- clusterProfiler::cnetplot(GeneList_Enrich_BP, foldChange = Data_filtered[GeneList_Matched$SYMBOL,2], circular = TRUE, colorEdge = TRUE), error=function(e) print("Not able to plot"))
    #Customized barplot
    tryCatch(GeneList_Enrich_BP_top <- as.data.frame(head(GeneList_Enrich_BP,ShowTermNum)), error=function(e) print("Not able to calculate"))
    tryCatch(GeneList_Enrich_BP_top$MinusLogPadj <- -log2(GeneList_Enrich_BP_top$p.adjust), error=function(e) print("Not able to calculate"))
    tryCatch(GeneList_Enrich_BP_top_sorted <- GeneList_Enrich_BP_top[order(GeneList_Enrich_BP_top$MinusLogPadj,decreasing = F),], error=function(e) print("Not able to calculate"))
    tryCatch(BP5 <- ggplot2::ggplot(data=GeneList_Enrich_BP_top_sorted,
                         ggplot2::aes(x=factor(Description,level=GeneList_Enrich_BP_top_sorted$Description),
                                      y=MinusLogPadj,fill=MinusLogPadj)) +
      ggplot2::scale_fill_gradient(low="blue", high="red", name = "- log2 (P adj)") +
      ggplot2::ggtitle(paste0("Gene ontology (GO) of ", OutputFileName)) +
      ggplot2::xlab("GO Terms, Biological Process (BP)\n\n") +
      ggplot2::ylab("\n\n- log2 (P adj)") +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 20, hjust = 0.5),
                     axis.text.x = ggplot2::element_text(size = 15),
                     axis.text.y = ggplot2::element_text(size = 15),
                     axis.title.x = ggplot2::element_text(size = 15),
                     axis.title.y = ggplot2::element_text(size = 15),
                     plot.margin = ggplot2::unit(c(1.5,1.5,1.5,1.5), units = "cm")) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::coord_flip(), error=function(e) print("Not able to plot"))
    tryCatch(BP6 <- clusterProfiler::heatplot(GeneList_Enrich_BP, foldChange=geneList, showCategory = ShowTermNum), error=function(e) print("Not able to plot"))
    pdf(paste0(OutputFileName,".GO.DEGenes.BP.pdf"),height = 10,width = 12)
    tryCatch(print(BP1), error=function(e) print("BP1 is missing"))
    tryCatch(print(BP2), error=function(e) print("BP2 is missing"))
    tryCatch(print(BP3), error=function(e) print("BP3 is missing"))
    tryCatch(print(BP4), error=function(e) print("BP4 is missing"))
    tryCatch(print(BP5), error=function(e) print("BP5 is missing"))
    tryCatch(print(BP6), error=function(e) print("BP6 is missing"))
    dev.off()

    print("Plotting figures for: CC")
    tryCatch(CC1 <- barplot(GeneList_Enrich_CC,showCategory = ShowTermNum), error=function(e) print("Not able to plot"))
    tryCatch(CC2 <- clusterProfiler::dotplot(GeneList_Enrich_CC,showCategory = ShowTermNum), error=function(e) print("Not able to plot"))
    tryCatch(CC3 <- clusterProfiler::cnetplot(GeneList_Enrich_CC), error=function(e) print("Not able to plot"))
    tryCatch(CC4 <- clusterProfiler::cnetplot(GeneList_Enrich_CC, foldChange = Data_filtered[GeneList_Matched$SYMBOL,2], circular = TRUE, colorEdge = TRUE), error=function(e) print("Not able to plot"))
    #Customized barplot
    tryCatch(GeneList_Enrich_CC_top <- as.data.frame(head(GeneList_Enrich_CC,ShowTermNum)), error=function(e) print("Not able to calculate"))
    tryCatch(GeneList_Enrich_CC_top$MinusLogPadj <- -log2(GeneList_Enrich_CC_top$p.adjust), error=function(e) print("Not able to calculate"))
    tryCatch(GeneList_Enrich_CC_top_sorted <- GeneList_Enrich_CC_top[order(GeneList_Enrich_CC_top$MinusLogPadj,decreasing = F),], error=function(e) print("Not able to calculate"))
    tryCatch(CC5 <- ggplot2::ggplot(data=GeneList_Enrich_CC_top_sorted,
                         ggplot2::aes(x=factor(Description,level=GeneList_Enrich_CC_top_sorted$Description),
                                      y=MinusLogPadj,fill=MinusLogPadj)) +
      ggplot2::scale_fill_gradient(low="blue", high="red", name = "- log2 (P adj)") +
      ggplot2::ggtitle(paste0("Gene ontology (GO) of ", OutputFileName)) +
      ggplot2::xlab("GO Terms, Cellular Component (CC)\n\n") +
      ggplot2::ylab("\n\n- log2 (P adj)") +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 20, hjust = 0.5),
                     axis.text.x = ggplot2::element_text(size = 15),
                     axis.text.y = ggplot2::element_text(size = 15),
                     axis.title.x = ggplot2::element_text(size = 15),
                     axis.title.y = ggplot2::element_text(size = 15),
                     plot.margin = ggplot2::unit(c(1.5,1.5,1.5,1.5), units = "cm")) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::coord_flip(), error=function(e) print("Not able to plot"))
    tryCatch(CC6 <- clusterProfiler::heatplot(GeneList_Enrich_CC, foldChange=geneList, showCategory = ShowTermNum), error=function(e) print("Not able to plot"))
    pdf(paste0(OutputFileName,".GO.DEGenes.CC.pdf"),height = 10,width = 12)
    tryCatch(print(CC1), error=function(e) print("CC1 is missing"))
    tryCatch(print(CC2), error=function(e) print("CC2 is missing"))
    tryCatch(print(CC3), error=function(e) print("CC3 is missing"))
    tryCatch(print(CC4), error=function(e) print("CC4 is missing"))
    tryCatch(print(CC5), error=function(e) print("CC5 is missing"))
    tryCatch(print(CC6), error=function(e) print("CC6 is missing"))
    dev.off()

    print("Plotting figures for: MF")
    tryCatch(MF1 <- barplot(GeneList_Enrich_MF,showCategory = ShowTermNum), error=function(e) print("Not able to plot"))
    tryCatch(MF2 <- clusterProfiler::dotplot(GeneList_Enrich_MF,showCategory = ShowTermNum), error=function(e) print("Not able to plot"))
    tryCatch(MF3 <- clusterProfiler::cnetplot(GeneList_Enrich_MF, foldChange = geneList), error=function(e) print("Not able to plot"))
    tryCatch(MF4 <- clusterProfiler::cnetplot(GeneList_Enrich_MF, foldChange = geneList, circular = TRUE, colorEdge = TRUE), error=function(e) print("Not able to plot"))
    #Customized barplot
    tryCatch(GeneList_Enrich_MF_top <- as.data.frame(head(GeneList_Enrich_MF,ShowTermNum)), error=function(e) print("Not able to calculate"))
    tryCatch(GeneList_Enrich_MF_top$MinusLogPadj <- -log2(GeneList_Enrich_MF_top$p.adjust), error=function(e) print("Not able to calculate"))
    tryCatch(GeneList_Enrich_MF_top_sorted <- GeneList_Enrich_MF_top[order(GeneList_Enrich_MF_top$MinusLogPadj,decreasing = F),], error=function(e) print("Not able to calculate"))
    tryCatch(MF5 <- ggplot2::ggplot(data=GeneList_Enrich_MF_top_sorted,
                         ggplot2::aes(x=factor(Description,level=GeneList_Enrich_MF_top_sorted$Description),
                                      y=MinusLogPadj,fill=MinusLogPadj)) +
      ggplot2::scale_fill_gradient(low="blue", high="red", name = "- log2 (P adj)") +
      ggplot2::ggtitle(paste0("Gene ontology (GO) of ", OutputFileName)) +
      ggplot2::xlab("GO Terms, Molecular Function (MF)\n\n") +
      ggplot2::ylab("\n\n- log2 (P adj)") +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 20, hjust = 0.5),
                     axis.text.x = ggplot2::element_text(size = 15),
                     axis.text.y = ggplot2::element_text(size = 15),
                     axis.title.x = ggplot2::element_text(size = 15),
                     axis.title.y = ggplot2::element_text(size = 15),
                     plot.margin = ggplot2::unit(c(1.5,1.5,1.5,1.5), units = "cm")) +
      ggplot2::geom_bar(stat="identity") +
      ggplot2::coord_flip(), error=function(e) print("Not able to plot"))
    tryCatch(MF6 <- clusterProfiler::heatplot(GeneList_Enrich_MF, foldChange=geneList, showCategory = ShowTermNum), error=function(e) print("Not able to plot"))
    pdf(paste0(OutputFileName,".GO.DEGenes.MF.pdf"),height = 10,width = 12)
    tryCatch(print(MF1), error=function(e) print("MF1 is missing"))
    tryCatch(print(MF2), error=function(e) print("MF2 is missing"))
    tryCatch(print(MF3), error=function(e) print("MF3 is missing"))
    tryCatch(print(MF4), error=function(e) print("MF4 is missing"))
    tryCatch(print(MF5), error=function(e) print("MF5 is missing"))
    tryCatch(print(MF6), error=function(e) print("MF6 is missing"))
    dev.off()
    print("Done GO analysis.")

    #GSEA
    print("Preparing a ranked gene list for GSEA")
    Data_log2FC <- Data$log2FoldChange
    names(Data_log2FC) <- Data$GeneName
    gene_list<-na.omit(Data_log2FC)
    gene_list = sort(gene_list, decreasing = TRUE)

    print("Starting GSEA analysis")
    gse <- gseGO(geneList=gene_list,
                 ont ="ALL",
                 keyType = "SYMBOL",
                 minGSSize = 3,
                 maxGSSize = 800,
                 pvalueCutoff = 1,
                 verbose = TRUE,
                 OrgDb = get(OrganDatabase),
                 pAdjustMethod = "BH")
    print("Writing GSEA results")
    write.table(gse, file = paste0(OutputFileName, ".GSEA.AllRankedGenes.txt"), sep = "\t", quote = F, col.names=NA)

    print("Plotting top 10 GSEA figures")
    pdf(paste0(OutputFileName, ".GSEA.Top10GSEA.pdf"))
    print("Please also check the GSEA.AllRankedGenes.txt file for their p values.")
    plot1 <- enrichplot::gseaplot2(gse, geneSetID = 1, title = gse@result$Description[1])
    plot2 <- enrichplot::gseaplot2(gse, geneSetID = 2, title = gse@result$Description[2])
    plot3 <- enrichplot::gseaplot2(gse, geneSetID = 3, title = gse@result$Description[3])
    plot4 <- enrichplot::gseaplot2(gse, geneSetID = 4, title = gse@result$Description[4])
    plot5 <- enrichplot::gseaplot2(gse, geneSetID = 5, title = gse@result$Description[5])
    plot6 <- enrichplot::gseaplot2(gse, geneSetID = 6, title = gse@result$Description[6])
    plot7 <- enrichplot::gseaplot2(gse, geneSetID = 7, title = gse@result$Description[7])
    plot8 <- enrichplot::gseaplot2(gse, geneSetID = 8, title = gse@result$Description[8])
    plot9 <- enrichplot::gseaplot2(gse, geneSetID = 9, title = gse@result$Description[9])
    plot10 <- enrichplot::gseaplot2(gse, geneSetID = 10, title = gse@result$Description[10])
    tryCatch(print(plot1), error=function(e) print("Not able to plot"))
    tryCatch(print(plot2), error=function(e) print("Not able to plot"))
    tryCatch(print(plot3), error=function(e) print("Not able to plot"))
    tryCatch(print(plot4), error=function(e) print("Not able to plot"))
    tryCatch(print(plot5), error=function(e) print("Not able to plot"))
    tryCatch(print(plot6), error=function(e) print("Not able to plot"))
    tryCatch(print(plot7), error=function(e) print("Not able to plot"))
    tryCatch(print(plot8), error=function(e) print("Not able to plot"))
    tryCatch(print(plot9), error=function(e) print("Not able to plot"))
    tryCatch(print(plot10), error=function(e) print("Not able to plot"))
    dev.off()
    print("Done GSEA analysis")
    rm(BP1,BP2,BP3,BP4,BP5,BP6)
    rm(CC1,CC2,CC3,CC4,CC5,CC6)
    rm(MF1,MF2,MF3,MF4,MF5,MF6)
    rm(GeneList_Enrich_BP_top_sorted,GeneList_Enrich_CC_top_sorted,GeneList_Enrich_MF_top_sorted)
    rm(plot1,plot2,plot3,plot4,plot5,plot6,plot7,plot8,plot9,plot10)

    print(paste0("Done running comparison: ",as.character(Compare_Sub$Condition_2)," v.s. ",as.character(Compare_Sub$Condition_1)))
    print("*****************************************************************")
    cat("\n")
  }
  print("Pipeline finished.")
}
