
##' EPIC_rel_F
##'
##' Function that generates the relative deconvolution plot with the EPIC method, with the fractional plot, and returns a df resulting from the deconvolution.
##' @param matrix expression matrix, this matrix must be in .txt format, organized in tabulated columns with one column containing the list of genes, with the name Gene.symbol, Gene.Symbol or Geneid, and in the rest of the columns the names of the samples and the expression values for the genes.  (Example. Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt).
##' @param sig.matrix Signature matrix, this matrix must be in .txt format, organized in tabulated columns with one column containing the list of genes, with the name Gene.symbol, Gene.Symbol or Geneid, and in the rest of the columns the cell types and the expression values for the genes.  (Example. Monaco-ABIS_Immune11_Microarray.txt).
##' @param fractions Vector that identifies the fractions of the sample list to be used to fractionate the samples in different graphs.
##' @param results_dir Directory to which the graphics will be generated.
##' @param height Number that defines the height of the graph, to increase or decrease the height, which has a default value of 10.
##' @param name Alternative name for the chart, which by default has the structure Sig_type_method.matrixName_Plot.png (Example. "graph_1")
##' @return Returns the generated graphs and the df of the deconvolution.
##' @author Nidia Barco Armengol
##' @export
##'
##' @examples
##' EPIC_rel_F("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", "Monaco-ABIS_Immune11_Microarray.txt", f, results_dir = "/bicoh/nidia/Deconv/Brain_MMartin/test", height= 8, name= "MONACO")

EPIC_rel_F <- function(matrix, sig.matrix, fractions, results_dir, height= 10, name= FALSE){
  require(EPIC)
  require(usethis)
  require(devtools)
  require(gplots)
  require(data.table)
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(immunedeconv)
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

  # Creació de la llista de gens i la matriu signature i generació de la matriu de deconvolució i conversió d'aquesta a df
  C_EPIC <- list()
  C_EPIC[["sigGenes"]] <- rownames(sig.mtrx)
  C_EPIC[["refProfiles"]] <- as.matrix(sig.mtrx)

  RESULTS <- t(EPIC::EPIC(bulk=cpm, reference=C_EPIC, withOtherCells=TRUE, scaleExprs=FALSE)$cellFractions) #scaleExprs=TRUE by default: only keep genes in common between matrices

  RESULTS <- as.data.frame(RESULTS)
  RESULTS$cell_type=rownames(RESULTS)
  E_rel <- RESULTS

  # Selecció de colors per als tipus celulars
  colors <- c(brewer.pal(12, "Set3"), brewer.pal(12, "Paired"))
  getPalette = colorRampPalette(colors)

  # Loop que genera una llista de llistes amb les posicions en que es troben les mostres de cada fracció
  list_pos <- list()
  unics <- unique(fractions)
  for (i in unics){
    positions <- c()
    for (x in 1:length(fractions)){
      element <- fractions[x]
      if (element==i){
        positions <- c(positions, x)
        list_pos[[i]] <- positions
      }
    }
  }
  # Creació del df de cada fracció i generació del seu grafic
  for (x in 1:length(list_pos)){
    e <- list_pos[x]
    df <- RESULTS$cell_type
    df <- as.data.frame(df)
    df_names <- c("cell_type")
    for (i in e){
      c <- RESULTS[,i]
      df <- cbind(df, c)
      df_names <- c(df_names, colnames(RESULTS)[i])
    }
    colnames(df) <- df_names

    if (name==FALSE){
      title <- gsub(".txt","_",sig.matrix)
      file <- c("EPIC_", title , x, "_","plot.png")
    }else{
      title <- name
      file <- c("EPIC_", title , "_", x, "_","plot.png")
    }
    file <- paste(file, collapse = "")

    p=df %>%
      gather(sample, fraction, -cell_type) %>%
      # plot as stacked bar chart
      ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
      geom_bar(stat='identity') +
      coord_flip() +
      scale_fill_manual(values = getPalette(length(rownames(RESULTS)))) +
      scale_x_discrete(limits = rev(levels(res_ciber)))
    ggsave(filename=paste(results_dir,file,sep="/"),dpi=300,width =9, height = height)
  }
  return(E_rel)
}
