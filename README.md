# easyRNAseqDE
An R package for RNAseq differential expression analysis and visualization, in an easy and reproducible way.

This package is for the downstream analysis after running [PipeRNAseq](https://github.com/sunyumail93/PipeRNAseq). Please download PipeRNAseq output files to the local computer and perform these analysis.

## Software prerequisites

This R package is designed to analyze RNAseq raw counts (produced by featureCounts from Subread) or TPM values from Salmon software.

The following packages in R are required

```
DESeq2
tximport
RColorBrewer
pheatmap
ggplot2
gplots
cluster
ggrepel
readr
reshape2
```

## Installation

```
library(devtools)
install_github("sunyumail93/easyRNAseqDE")

library(easyRNAseqDE)

data(featureCounts_count_matrix)

head(featureCounts_count_matrix)

# Output:
              neg1 neg2 neg3 neg4 pos_1 pos_2 pos_3 control_1 control_2
RP23-271O17.1    0    0    0    0     0     0     0         0         0
Gm26206          0    0    0    0     0     0     0         0         0
Xkr4             0    0    0    0     6     3     1         0         0
RP23-317L18.1    0    0    0    0     0     0     0         0         0
RP23-317L18.4    0    0    0    0     0     0     0         0         0
RP23-317L18.3    0    0    0    0     0     0     0         0         0
```


## 1, Merge featureCounts results into a table/matrix

This function is compatible with the featureCounts output files from PipeRNAseq, with *.featureCounts.gene.txt suffix.

```
# Prepare the featureCounts output files after running PipeRNAseq
Control.rep1.featureCounts.gene.txt
Control.rep2.featureCounts.gene.txt
Treated.rep1.featureCounts.gene.txt
Treated.rep2.featureCounts.gene.txt

# Prepare the sample information file, with contains 3 columns: SampleFile, ShortName, Condition
# The SampleFile should be exactly same as the file names above
# The ShortName indicates the name in the final merged table, to make the names simplified
# The Condition column is for future DE analysis, but won't be used in this function

$ cat SampleInfoFile.txt
Control.rep1.featureCounts.gene.txt Control.rep1 Control
Control.rep2.featureCounts.gene.txt Control.rep2 Control
Treated.rep1.featureCounts.gene.txt Treated.rep1 Treated
Treated.rep2.featureCounts.gene.txt Treated.rep2 Treated

# Then we can run the following function:
# Please use setwd() to navigate to the working directory containing featureCounts output and SampleInfoFile.txt
MergeFeatureCounts(SampleInfoFile, OutputPrefix="default", write_to_file=F)

# Example 1: return a merged table as data.frame
MergedCounts <- MergeFeatureCounts(SampleInfoFile = "SampleInfoFile.txt")

# Example 2, return a merged table as data.frame, and write to a file
MergedCounts <- MergeFeatureCounts(SampleInfoFile = "SampleInfoFile.txt", 
                   OutputPrefix = "Merged.featureCounts.txt",
                   write_to_file = F)
```

This MergedCounts table can be directly used in the DE analysis, as the input data for **easyDE_FromRawCounts**.

## 2, Merge salmon results into a table/matrix

This function is compatible with the salmon output files from PipeRNAseq, with *.salmon.sf suffix.

```
# Prepare the featureCounts output files after running PipeRNAseq
Control.rep1.quant.sf
Control.rep2.quant.sf
Treated.rep1.quant.sf
Treated.rep2.quant.sf

# Copy the uniqueMatchingFile file from PipeRNAseq pipeline
mm10.uniqueMatching.txt  #For mouse
hg38.uniqueMatching.txt  #For human

# Prepare the sample information file, with contains 3 columns: SampleFile, ShortName, Condition
# The SampleFile should be exactly same as the file names above
# The ShortName indicates the name in the final merged table, to make the names simplified
# The Condition column is for future DE analysis, but won't be used in this function

$ cat SampleInfoFile.txt
Control.rep1.quant.sf Control.rep1 Control
Control.rep2.quant.sf Control.rep2 Control
Treated.rep1.quant.sf Treated.rep1 Treated
Treated.rep2.quant.sf Treated.rep2 Treated

# Then we can run the following function:
# Please use setwd() to navigate to the working directory containing featureCounts output and SampleInfoFile.txt
MergeSalmon(SampleInfoFile, uniqueMatchingFile, OutputPrefix="default", write_to_file=F)

# Example 1: return a merged Counts and TPM table as data.frame
MergedCountsTPM <- MergeSalmon(SampleInfoFile = "SampleInfoFile.txt", 
                   uniqueMatchingFile = "mm10.uniqueMatching.txt")

# Example 2, return a merged Counts and TPM table as data.frame, and write to a file
MergedCountsTPM <- MergeSalmon(SampleInfoFile = "SampleInfoFile.txt", 
                   uniqueMatchingFile = "mm10.uniqueMatching.txt",
                   OutputPrefix = "Merged.salmon.CountsAndTPM.txt",
                   write_to_file = T)
```

This merged salmon table can be used to examine the TPM values, but won't be compatible with the DE analysis. To perform DE analysis using salmon outputs, please refer section 4 and use **easyDE_FromSalmon** function.

## 3, Run differential expression (DE) analysis from merged featureCounts table

We can use any merged raw counts table (either from **MergeFeatureCounts**, or read from an external file) and perform DE analysis using one of the main function in this package: **easyDE_FromRawCounts**

```
# Prepare the sample information/Label file, with contains 3 columns: SampleFile, ShortName, Condition
# The SampleFile should be exactly same as the file names above
# The ShortName indicates the name in the final merged table, to make the names simplified
# The Condition column is for future DE analysis, but won't be used in this function

$ cat SampleInfoFile.txt
PS.rep1.featureCounts.gene.txt PS.rep1 PS
PS.rep2.featureCounts.gene.txt PS.rep2 PS
RS.rep1.featureCounts.gene.txt RS.rep1 RS
RS.rep2.featureCounts.gene.txt RS.rep2 RS
ART.rep1.featureCounts.gene.txt ART.rep1 ART
ART.rep2.featureCounts.gene.txt ART.rep2 ART

# Prepare the comparison file, with 3 columns: Condition1, Condition2, ComparisonName/Output Prefix
# The direction of the comparison: Condition2_vs_Condition1:

$ cat comparison.txt
PS RS Output.Human.RSvsPS
PS ART Output.Human.ARTvsPS

# Prepare the raw counts using MergeFeatureCounts function mentioned above
MergedCounts <- MergeFeatureCounts(SampleInfoFile = "SampleInfoFile.txt")

# Turn on the createQuickomicsFiles=T and provide a output prefix for Quickomics file generation

# Run the full DE analysis
easyDE_FromRawCounts(count_matrix = MergedCounts, 
                     LabelFile = "SampleInfoFile.txt", 
                     ComparisonFile = "comparison.txt", 
                     createQuickomicsFiles=T, 
                     QuickomicsPrefix="Two_DE_Analysis")
```

For each comparison, the following files will be generated. Here, the Prefix = "Output.Human.RSvsPS" or "Output.Human.ARTvsPS":

```
Prefix.DESeq2.txt
Prefix.DistancePlot.pdf
Prefix.HeatmapTop50.pdf
Prefix.MAPlot.Labeltop20.pdf
Prefix.PCA.1.pdf
Prefix.PCA.2.pdf
Prefix.VolcanoPlot.pdf
```

[Quickomics](https://github.com/interactivereport/Quickomics) files for the whole project:

```
Two_DE_Analysis_Comparison_data.csv
Two_DE_Analysis_Exp_Log2Data.csv             #Choose one from Log2Data and RawData to upload to Quickomics
Two_DE_Analysis_Exp_RawData.csv
Two_DE_Analysis_ProteinGeneName_optional.csv
Two_DE_Analysis_Sample_metadata.csv
```

## 4, Run differential expression (DE) analysis from salmon output

This function is compatible with the salmon output files from PipeRNAseq, with *.salmon.sf suffix. Different from **MergeSalmon** function, **easyDE_FromSalmon** runs the full DE analysis.

```
# Prepare the featureCounts output files after running PipeRNAseq
Control.rep1.quant.sf
Control.rep2.quant.sf
Treated.rep1.quant.sf
Treated.rep2.quant.sf

# Copy the uniqueMatchingFile file from PipeRNAseq pipeline
mm10.uniqueMatching.txt  #For mouse
hg38.uniqueMatching.txt  #For human

# Prepare the sample information file, with contains 3 columns: SampleFile, ShortName, Condition
# The SampleFile should be exactly same as the file names above
# The ShortName indicates the name in the final merged table, to make the names simplified
# The Condition column is for future DE analysis, but won't be used in this function

$ cat SampleInfoFile.txt
Control.rep1.quant.sf Control.rep1 Control
Control.rep2.quant.sf Control.rep2 Control
Treated.rep1.quant.sf Treated.rep1 Treated
Treated.rep2.quant.sf Treated.rep2 Treated

# Prepare the comparison file, with 3 columns: Condition1, Condition2, ComparisonName/Output Prefix
# The direction of the comparison: Condition2_vs_Condition1:
# You can include multiple lines of comparisons here. For each comparison, no need to cover all data.
# The pipeline will subset data included in the comparison and run the analysis with those data only.

$ cat comparison.txt
Control Treated Output.Control_vs_Treated

# Then we can run the following function:
# Please use setwd() to navigate to the working directory containing featureCounts output and SampleInfoFile.txt

# Example 1: return a merged Counts and TPM table as data.frame
easyDE_FromSalmon(SampleInfo = "SampleInfoFile.txt",
                  uniqueMatchingFile = "mm10.uniqueMatching.txt",
                  ComparisonFile = "comparison.txt", 
                  createQuickomicsFiles=T, QuickomicsPrefix="DEAnalysis")
```

For each comparison, the following files will be generated. Here, the Prefix = "Output.Control_vs_Treated":

```
Prefix.salmon.All.DEseq2GeneCountsTPMrlogNormcounts.txt
Prefix.salmon.Plot.DistancePlot.All.pdf
Prefix.salmon.Plot.HeatmapTop50.pdf
Prefix.salmon.Plot.MA.pdf
Prefix.salmon.Plot.PCA.1.pdf
Prefix.salmon.Plot.PCA.2.pdf
Prefix.salmon.Plot.VolcanoPlot.pdf
```

[Quickomics](https://github.com/interactivereport/Quickomics) files for the whole project:

```
DEAnalysis_Comparison_data.csv
DEAnalysis_Exp_Log2TPMData.csv             #Choose one from Log2Data and RawData to upload to Quickomics
                                           #log2(TPM+1) is good for plotting heatmap/PCA plots, without any other adjustments.
DEAnalysis_Exp_RawTPMData.csv              #Raw TPM is suggested if you woule like to plot the expression values.
DEAnalysis_Exp_rlogData.csv                #rlog is better at plotting heatmap/PCA plots, and it is adjusted for variance stablization
DEAnalysis_ProteinGeneName_optional.csv
DEAnalysis_Sample_metadata.csv
```
