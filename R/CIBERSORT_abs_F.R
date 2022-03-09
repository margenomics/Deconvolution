
##' CIBERSORT_abs_F
##'
##' Function that generates the absolute deconvolution plot with the CIBERSORT method, with the fractional plot, and returns a df resulting from the deconvolution.
##' @param matrix expression matrix, this matrix must be in .txt format, organized in tabulated columns with one column containing the list of genes, with the name Gene.symbol, Gene.Symbol or Geneid, and in the rest of the columns the names of the samples and the expression values for the genes.  (Example. Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt).
##' @param sig.matrix Signature matrix, this matrix must be in .txt format, organized in tabulated columns with one column containing the list of genes, with the name Gene.symbol, Gene.Symbol or Geneid, and in the rest of the columns the cell types and the expression values for the genes.  (Example. Monaco-ABIS_Immune11_Microarray.txt).
##' @param fractions Vector that identifies the fractions of the sample list to be used to fractionate the samples in different graphs.
##' @param height Number that defines the height of the graph, to increase or decrease the height, which has a default value of 10.
##' @param results_dir Directory to which the graphics will be generated.
##' @param name Alternative name for the chart, which by default has the structure Sig_type_method.matrixName_Plot.png (Example. "graph_1")
##' @return Returns the generated graphs and the df of the deconvolution.
##' @author Nidia Barco Armengol
##' @export
##'
##' @examples
##' CIBERSORT_abs_F("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", "Monaco-ABIS_Immune11_Microarray.txt", f, results_dir = "/bicoh/nidia/Deconv/Brain_MMartin/test", height= 8, name, "MONACO")

CIBERSORT_abs_F <- function(matrix, sig.matrix, fractions, results_dir, height= 10, name= FALSE){
  require(usethis)
  require(devtools)
  require(gplots)
  require(data.table)
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  require(immunedeconv)
  require(RColorBrewer)
  source("CIBERSORT.R")
  # Construcció de la matriu de deconvolució absoluta i conversió d'aquesta en un dataframe
  res_ciber = CIBERSORT(sig_matrix=sig.matrix,
                        mixture_file = matrix,
                        QN = F,perm=100,absolute = T,abs_method='no.sumto1')
  res_ciber = as.data.frame(t(res_ciber))
  res_ciber$cell_type=rownames(res_ciber)
  res_ciber$cell_type <- NULL
  C_abs <- res_ciber[1:(length(rownames(res_ciber))-4),]

  # Selecció de colors ajuntant dues paletes de colors de colorbrewer
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
    df <- rownames(res_ciber)
    df <- as.data.frame(df)
    df_names <- c("cell_type")
    for (i in e){
      c <- res_ciber[,i]
      df <- cbind(df, c)
      df_names <- c(df_names, colnames(res_ciber)[i])
    }
    colnames(df) <- df_names

    if (name==FALSE){
      title <- gsub(".txt","_",sig.matrix)
      file <- c("CIBERSORT.abs_", title , x, "_","plot.png")
    }else{
      title <- name
      file <- c("CIBERSORT.abs_", title , "_", x, "_","plot.png")
    }
    file <- paste(file, collapse = "")

    p=df[1:((grep("P-value",rownames(res_ciber)))-1),] %>%
      gather(sample, fraction, -cell_type) %>%
      # plot as stacked bar chart
      ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
      geom_bar(stat='identity') +
      coord_flip() +
      scale_fill_manual(values = getPalette(length(rownames(res_ciber)))) +
      scale_x_discrete(limits = rev(levels(res_ciber)))
    ggsave(filename=paste(results_dir,file,sep="/"),dpi=300,width =9, height = height)
  }
  return(C_abs)
}
