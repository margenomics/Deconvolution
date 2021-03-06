% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Deconv.R
\name{deconv}
\alias{deconv}
\title{DECONVOLUTION function}
\usage{
deconv(
  matrix,
  sig.matrix,
  method,
  type,
  fraccionate.samples = FALSE,
  fractions = NULL,
  results_dir,
  height = NULL,
  name = FALSE
)
}
\arguments{
\item{matrix}{expression matrix, this matrix must be in .txt format, organized in tabulated columns with one column containing the list of genes, with the name Gene.symbol, Gene.Symbol or Geneid, and in the rest of the columns the names of the samples and the expression values for the genes.  (Example. Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt).}

\item{sig.matrix}{Signature matrix, this matrix must be in .txt format, organized in tabulated columns with one column containing the list of genes, with the name Gene.symbol, Gene.Symbol or Geneid, and in the rest of the columns the cell types and the expression values for the genes.  (Example. Monaco-ABIS_Immune11_Microarray.txt).}

\item{method}{keyword to choose the deconvolution method, accept "CIBERSORT", "EPIC", "FARDEEP" and "ALL".}

\item{type}{Keyword allowing to choose the type of evolution, relative ("rel"), absolute ("abs") or both ("abs_rel").}

\item{fraccionate.samples}{Parameter that accepts TRUE or FALSE to split or not the samples of the expression matrix in different graphs.}

\item{fractions}{Vector that identifies the fractions of the sample list to be used to fractionate the samples in different graphs.}

\item{results_dir}{Directory to which the graphics will be generated.}

\item{height}{Number that defines the height of the graph, to increase or decrease the height, which has a default value of 10.}

\item{name}{Alternative name for the chart, which by default has the structure Sig_type_method.matrixName_Plot.png (Example. "graph_1")}
}
\value{
The function returns the deconvolution graphics in .png format to the working directory, and a list of the df's resulting from the deconvolution.
}
\description{
Function that generates massive deconvolution plots with the CIBERSORT, EPIC and FARDEEP methods, and returns a list of df's resulting from the deconvolutions.
}
\examples{
deconv("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", "Monaco-ABIS_Immune11_Microarray.txt", "ALL", "rel", fraccionate.samples=TRUE, fractions= f, results_dir = "/bicoh/nidia/Deconv/Brain_MMartin/test", height= 8, name= "MONACO")
}
\author{
Nidia Barco Armengol
}
