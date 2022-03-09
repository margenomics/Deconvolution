
##' FARDEEP_abs_NF
##'
##' Function that generates the absolute deconvolution plot with the FARDEEP method, with the non-fractional plot, and returns a df resulting from the deconvolution.
##' @param matrix expression matrix, this matrix must be in .txt format, organized in tabulated columns with one column containing the list of genes, with the name Gene.symbol, Gene.Symbol or Geneid, and in the rest of the columns the names of the samples and the expression values for the genes.  (Example. Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt).
##' @param sig.matrix Signature matrix, this matrix must be in .txt format, organized in tabulated columns with one column containing the list of genes, with the name Gene.symbol, Gene.Symbol or Geneid, and in the rest of the columns the cell types and the expression values for the genes.  (Example. Monaco-ABIS_Immune11_Microarray.txt).
##' @param results_dir Directory to which the graphics will be generated.
##' @param height Number that defines the height of the graph, to increase or decrease the height, which has a default value of 10.
##' @param name Alternative name for the chart, which by default has the structure Sig_type_method.matrixName_Plot.png (Example. "graph_1")
##' @return Returns the generated graph and the df of the deconvolution.
##' @author Nidia Barco Armengol
##' @export
##'
##' @examples
##' FARDEEP_abs_NF("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", "Monaco-ABIS_Immune11_Microarray.txt", results_dir = "/bicoh/nidia/Deconv/Brain_MMartin/test", height= 8, name= "MONACO")

FARDEEP_abs_NF <- function(matrix, sig.matrix, results_dir, height= 10, name= FALSE){
  require(usethis)
  require(devtools)
  require(gplots)
  require(data.table)
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(immunedeconv)
  require(FARDEEP)
  require(RColorBrewer)
  # Llegeix la matriu d'expressió en txt i posa els rownames amb els noms de la columna de gens, 3 possibles noms
  cpm=read.table(matrix,sep = "\t",header = T)

  if (colnames(cpm)[1]=="Gene.symbol"){
    cpm <- cpm %>%
      group_by(Gene.symbol) %>%
      summarise_all(mean)
    names <- cpm$Gene.symbol
    cpm=as.matrix(cpm[,-1])
    rownames(cpm)=names
  }
  if (colnames(cpm)[1]=="Gene.Symbol"){
    cpm <- cpm %>%
      group_by(Gene.Symbol) %>%
      summarise_all(mean)
    names <- cpm$Gene.Symbol
    cpm=as.matrix(cpm[,-1])
    rownames(cpm)=names
  }
  if (colnames(cpm)[1]=="Geneid"){
    cpm <- cpm %>%
      group_by(Geneid) %>%
      summarise_all(mean)
    names <- cpm$Geneid
    cpm=as.matrix(cpm[,-1])
    rownames(cpm)=names
  }

  # Llegeix la matriu signature en txt i posa els rownames amb els noms de la columna de gens, 3 possibles noms
  sig.mtrx=read.table(sig.matrix,sep = "\t",header = T)

  if (colnames(sig.mtrx)[1]=="Gene.symbol"){
    sig.mtrx <- sig.mtrx %>%
      group_by(Gene.symbol) %>%
      summarise_all(mean)
    names <- sig.mtrx$Gene.symbol
    sig.mtrx=as.matrix(sig.mtrx[,-1])
    rownames(sig.mtrx)=names
  }
  if (colnames(sig.mtrx)[1]=="Gene.Symbol"){
    sig.mtrx <- sig.mtrx %>%
      group_by(Gene.Symbol) %>%
      summarise_all(mean)
    names <- sig.mtrx$Gene.Symbol
    sig.mtrx=as.matrix(sig.mtrx[,-1])
    rownames(sig.mtrx)=names
  }
  if (colnames(sig.mtrx)[1]=="Geneid"){
    sig.mtrx <- sig.mtrx %>%
      group_by(Geneid) %>%
      summarise_all(mean)
    names <- sig.mtrx$Geneid
    sig.mtrx=as.matrix(sig.mtrx[,-1])
    rownames(sig.mtrx)=names
  }

  # Construcció de la matriu de deconvolució absoluta i conversió d'aquesta en un dataframe
  RESULTS = t(FARDEEP::fardeep(sig.mtrx, cpm, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)$abs.beta)
  #RESULTS = t(FARDEEP::fardeep(sig.mtrx, cpm, nn = TRUE, intercept = TRUE, permn = 10, QN = FALSE)$relative.beta)

  RESULTS <- as.data.frame(RESULTS)
  RESULTS$cell_type=rownames(RESULTS)
  F_abs <- RESULTS

  # Defineix el nom de l'arxiu depenent de si l'usuari introdueix el nom o si es predeterminat
  if (name==FALSE){
    title <- gsub(".txt","_",sig.matrix)
    file <- c("FARDEEP.abs_", title , "_","plot.png")
  }else{
    title <- name
    file <- c("FARDEEP.abs_", title , "_","plot.png")
  }
  file <- paste(file, collapse = "")

  # Selecció de colors per als tipus celulars
  colors <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
  getPalette = colorRampPalette(colors)

  # Creació del grafic i guardat d'aquest al directori de resultats
  p=RESULTS %>%
    gather(sample, fraction, -cell_type) %>%
    # plot as stacked bar chart
    ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
    geom_bar(stat='identity') +
    coord_flip() +
    scale_fill_manual(values = getPalette(length(rownames(RESULTS)))) +
    scale_x_discrete(limits = rev(levels(RESULTS)))
  ggsave(filename=paste(results_dir,file,sep="/"),dpi=300,width =9, height = height)
  return(F_abs)
}
