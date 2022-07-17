#' Run KEGG analysis using DE results.
#' @param ComparisonFile DE comparison file with three columns: Condition_1 Condition_2 OutputPrefix.
#' @param OrganDatabase Organism GO annotation database, for example, 'org.Mm.eg.db' for mouse, and org.Hs.eg.db for human.
#' @param KEGGSpeciesName KEGG pathway species name, for example, 'mmu' for mouse, 'hsa' for human.
#' @param KEGGPathwayFile A two-column file for KEGG pathways to plot: KEGG_id KEGG_full_name, default="NA"
#' @param GOpajd_cutoff pajd cutoff to define significant genes, default=0.05.
#' @param GOlog2FC_cutoff log2 fold change cutoff to define significant genes, default=0.05
#' @param forcerun Force to run the analysis if DE gene number<10, default=False.
#' @export
#
# Run KEGG analysis using DE results. Please run this function after running easyDE_FromRawCounts or easyDE_FromSalmon.
# The output results include KEGG enrichment (using significant DE genes) and GSEA results (using all ranked genes).
# The output also includes a set of KEGG figures based on given KEGG pathway IDs.

easyEnrich_KEGG <- function(ComparisonFile, OrganDatabase, KEGGSpeciesName, KEGGPathwayFile=NA, GOpajd_cutoff=0.05, GOlog2FC_cutoff=1.5, forcerun=F){
  print("Start running easyGO_KEGG")
  Comparison <- read.table(ComparisonFile,header=F,col.names = c("Condition_1","Condition_2","OutputFileName"))

  if (is.na(KEGGPathwayFile)){
    print("No KEGG pathway file provided. No figures will be plotted.")
  }else{
    print("Reading KEGG file...")
    KEGG_id <- read.table(KEGGPathwayFile,header=F,col.names = c("PathwayID","Name"))
  }

  for (i in 1:dim(Comparison)[1]){
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

    #KEGG enrichment analysis
    print("Running KEGG enrichment analysis using significant DE genes...")
    kk <- clusterProfiler::enrichKEGG(gene         = GeneList_Matched$ENTREZID,
                                      organism     = KEGGSpeciesName,
                                      pvalueCutoff = 1)
    write.table(kk, file = paste0(OutputFileName, ".KEGG.DEGenes.txt"), sep = "\t", quote = F, col.names=NA)

    print("Preparing a ranked gene list for GSEA")
    Data_log2FC <- Data$log2FoldChange
    names(Data_log2FC) <- Data$GeneName
    gene_list<-na.omit(Data_log2FC)
    gene_list = sort(gene_list, decreasing = TRUE)
    GeneList_Matched = clusterProfiler::bitr(as.character(names(gene_list)), fromType="SYMBOL", toType="ENTREZID", OrgDb=OrganDatabase)
    gene_list_filtered <- gene_list[GeneList_Matched$SYMBOL]
    names(gene_list_filtered) <- GeneList_Matched$ENTREZID
    print("Starting KEGG GSEA analysis using ranked gene list")
    kk2 <- clusterProfiler::gseKEGG(geneList     = gene_list_filtered,
                                    organism     = 'mmu',
                                    pvalueCutoff = 1,
                                    verbose      = F)
    write.table(kk2, file = paste0(OutputFileName, ".KEGG.AllRankedGenes.txt"), sep = "\t", quote = F, col.names=NA)

    if (!is.na(KEGGPathwayFile)){
      print("Start plotting pathway maps...")
      dir.create(paste0(OutputFileName, ".KEGG.Maps"))
      wd <- getwd()
      setwd(paste0(getwd(),"/",OutputFileName, ".KEGG.Maps"))
      for (j in 1:dim(KEGG_id)[1]){
        print(j)
        CurrentPathway=KEGG_id[j,]$PathwayID
        CurrentPathwayName=KEGG_id[j,]$Name
        print(paste0("Generating pathway map ", j," for: ", CurrentPathway, ", ", CurrentPathwayName))
        detach("package:pathfindR", unload=TRUE)
        detach("package:pathview", unload=TRUE)
        pathview::pathview(gene.data  = gene_list_filtered,
                           pathway.id = CurrentPathway,
                           species    = KEGGSpeciesName,
                           out.suffix=OutputFileName)
      }
      setwd(wd)
    }
    print("*****************************************************************")
    cat("\n")
  }
  print("Pipeline finished.")
}
