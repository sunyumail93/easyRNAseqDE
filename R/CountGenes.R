#' Count significantly changed genes from DESeq2 output
#'
#' @param res_ordered sorted DESeq2 result table.
#' @return A summary of significantly changed genes, based on padj<0.05 and padj<0.01 cutoffs.
#' @export
#
# Count significantly changed genes

CountGenes <- function(res_ordered){
  #Define significant genes
  print("Counting significanly changed genes...")
  res_ordered_005 <- subset(res_ordered,res_ordered$padj<0.05)
  res_ordered_005Up <- res_ordered_005[res_ordered_005$log2FoldChange > 0,]
  res_ordered_005Down <- res_ordered_005[res_ordered_005$log2FoldChange < 0,]
  print(paste(dim(res_ordered_005)[1]," significantly changed genes with padj < 0.05 cutoff, ",dim(res_ordered_005Up)[1]," up-regulated, ",dim(res_ordered_005Down)[1]," down-regulated",sep=""))

  res_ordered_001 <- subset(res_ordered,res_ordered$padj<0.01)
  res_ordered_001Up <- res_ordered_001[res_ordered_001$log2FoldChange > 0,]
  res_ordered_001Down <- res_ordered_001[res_ordered_001$log2FoldChange < 0,]
  print(paste(dim(res_ordered_001)[1]," Genes with padj < 0.01 cutoff, ",dim(res_ordered_001Up)[1]," up-regulated, ",dim(res_ordered_001Down)[1]," down-regulated",sep=""))
}
