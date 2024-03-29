% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/easyEnrich_KEGG.R
\name{easyEnrich_KEGG}
\alias{easyEnrich_KEGG}
\title{Run KEGG analysis using DE results.}
\usage{
easyEnrich_KEGG(
  ComparisonFile,
  OrganDatabase,
  KEGGSpeciesName,
  KEGGPathwayFile = NA,
  GOpajd_cutoff = 0.05,
  GOlog2FC_cutoff = 1.5,
  forcerun = F
)
}
\arguments{
\item{ComparisonFile}{DE comparison file with three columns: Condition_1 Condition_2 OutputPrefix.}

\item{OrganDatabase}{Organism GO annotation database, for example, 'org.Mm.eg.db' for mouse, and org.Hs.eg.db for human.}

\item{KEGGSpeciesName}{KEGG pathway species name, for example, 'mmu' for mouse, 'hsa' for human.}

\item{KEGGPathwayFile}{A two-column file for KEGG pathways to plot: KEGG_id KEGG_full_name, default="NA"}

\item{GOpajd_cutoff}{pajd cutoff to define significant genes, default=0.05.}

\item{GOlog2FC_cutoff}{log2 fold change cutoff to define significant genes, default=1.5.}

\item{forcerun}{Force to run the analysis if DE gene number<10, default=False.}
}
\value{
For each comparison, a KEGG summary table, and a folder containing KEGG plots using pathway names provided by KEGGPathwayFile will be generated.
}
\description{
Run KEGG analysis using DE results.
}
