
##' DECONVOLUTION function
##'
##' Function that generates massive deconvolution plots with the CIBERSORT, EPIC and FARDEEP methods, and returns a list of df's resulting from the deconvolutions.
##' @param matrix expression matrix, this matrix must be in .txt format, organized in tabulated columns with one column containing the list of genes, with the name Gene.symbol, Gene.Symbol or Geneid, and in the rest of the columns the names of the samples and the expression values for the genes.  (Example. Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt).
##' @param sig.matrix Signature matrix, this matrix must be in .txt format, organized in tabulated columns with one column containing the list of genes, with the name Gene.symbol, Gene.Symbol or Geneid, and in the rest of the columns the cell types and the expression values for the genes.  (Example. Monaco-ABIS_Immune11_Microarray.txt).
##' @param method keyword to choose the deconvolution method, accept "CIBERSORT", "EPIC", "FARDEEP" and "ALL".
##' @param type Keyword allowing to choose the type of evolution, relative ("rel"), absolute ("abs") or both ("abs_rel").
##' @param fraccionate.samples Parameter that accepts TRUE or FALSE to split or not the samples of the expression matrix in different graphs.
##' @param fractions Vector that identifies the fractions of the sample list to be used to fractionate the samples in different graphs.
##' @param results_dir Directory to which the graphics will be generated.
##' @param height Number that defines the height of the graph, to increase or decrease the height, which has a default value of 10.
##' @param name Alternative name for the chart, which by default has the structure Sig_type_method.matrixName_Plot.png (Example. "graph_1")
##' @return The function returns the deconvolution graphics in .png format to the working directory, and a list of the df's resulting from the deconvolution.
##' @author Nidia Barco Armengol
##' @export
##'
##' @examples
##' deconv("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", "Monaco-ABIS_Immune11_Microarray.txt", "ALL", "rel", fraccionate.samples=TRUE, fractions= f, results_dir = "/bicoh/nidia/Deconv/Brain_MMartin/test", height= 8, name= "MONACO")


deconv <- function(matrix, sig.matrix, method, type, fraccionate.samples= FALSE, fractions= NULL, results_dir, height= NULL, name= FALSE){
  require(EPIC)
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
  source("CIBERSORT.R")
  # Primera selecció, entra si NO es volen generar grafics fraccionats
  if (fraccionate.samples==FALSE){
    # Filtre per fer tots els metodes
    if (method=="ALL"){
      # Filtre per fer només el gràfic de deconvolució relativa
      if (type=="rel"){
        C_rel <- CIBERSORT_rel_NF(matrix, sig.matrix, results_dir, height, name)
      }
      # Filtre per fer només el grafic de deconvolució absoluta
      if (type=="abs"){
        C_abs <- CIBERSORT_abs_NF(matrix, sig.matrix, results_dir, height, name)
      }
      # Filtre per fer els grafics de deconvolució relativa i absoluta
      if (type=="abs_rel"){
        C_rel <- CIBERSORT_rel_NF(matrix, sig.matrix, results_dir, height, name)
        C_abs <- CIBERSORT_abs_NF(matrix, sig.matrix, results_dir, height, name)
      }
      # Filtre per fer només el gràfic de deconvolució relativa
      if (type=="rel"){
        E_rel <- EPIC_rel_NF(matrix, sig.matrix, results_dir, height, name)
      }
      # Filtre per fer només el gràfic de deconvolució absoluta, en aquest cas com EPIC no te deconvolució absoluta mostra un warning informatiu
      if (type=="abs"){
        warning("Warning: EPIC only has relative type.")
      }
      # Filtre per fer els grafics de deconvolució relativa i absoluta, mostra el warning per absoluta
      if (type=="abs_rel"){
        E_rel <- EPIC_rel_NF(matrix, sig.matrix, results_dir, height, name)
        warning("Warning: EPIC only has relative type.")
      }
      # Filtre per fer només el gràfic de deconvolució relativa
      if (type=="rel"){
        F_rel <- FARDEEP_rel_NF(matrix, sig.matrix, results_dir, height, name)
        df_list <- list(C_rel, E_rel, F_rel)
        names(df_list) <- c("C_rel", "E_rel", "F_rel")
      }
      # Filtre per fer només el gràfic de deconvolució absoluta
      if (type=="abs"){
        F_abs <- FARDEEP_abs_NF(matrix, sig.matrix, results_dir, height, name)
        df_list <- list(C_abs, F_abs)
        names(df_list) <- c("C_abs", "F_abs")
      }
      # Filtre per fer els grafics de deconvolució relativa i absoluta
      if (type=="abs_rel"){
        F_rel <- FARDEEP_rel_NF(matrix, sig.matrix, results_dir, height, name)
        F_abs <- FARDEEP_abs_NF(matrix, sig.matrix, results_dir, height, name)
        df_list <-list(C_rel, C_abs, E_rel, F_rel, F_abs)
        names(df_list) <- c("C_rel", "C_abs", "E_rel", "F_rel", "F_abs")
      }
    }
    # Filtre per generar només els grafics de CIBERSORT
    if (method=="CIBERSORT"){
      # Filtre per fer només el gràfic de deconvolució relativa
      if (type=="rel"){
        C_rel <- CIBERSORT_rel_NF(matrix, sig.matrix, results_dir, height, name)
        df_list <- list(C_rel)
        names(df_list) <- c("C_rel")
      }
      # Filtre per fer només el gràfic de deconvolució absoluta
      if (type=="abs"){
        C_abs <- CIBERSORT_abs_NF(matrix, sig.matrix, results_dir, height, name)
        df_list <- list(C_abs)
        names(df_list) <- c("C_abs")
      }
      # Filtre per fer els grafics de deconvolució relativa i absoluta
      if (type=="abs_rel"){
        C_rel <- CIBERSORT_rel_NF(matrix, sig.matrix, results_dir, height, name)
        C_abs <- CIBERSORT_abs_NF(matrix, sig.matrix, results_dir, height, name)
        df_list <- list(C_rel, C_abs)
        names(df_list) <- c("C_rel", "C_abs")
      }
    }
    # Filtre per generar només els grafics de EPIC
    if (method=="EPIC"){
      # Filtre per fer només el gràfic de deconvolució relativa
      if (type=="rel"){
        E_rel <- EPIC_rel_NF(matrix, sig.matrix, results_dir, height, name)
        df_list <-list(E_rel)
        names(df_list) <- c("E_rel")
      }
      # Filtre per fer només el gràfic de deconvolució absoluta, en aquest cas posa warning informatiu
      if ((type=="abs")){
        warning("Warning: EPIC only has relative type.")
        df_list <- list()
      }
      # Filtre per fer els grafics de deconvolució relativa i absoluta, per absoluta un warning
      if (type=="abs_rel"){
        E_rel <- EPIC_rel_NF(matrix, sig.matrix, results_dir, height, name)
        warning("Warning: EPIC only has relative type.")
        df_list <- list(E_rel)
        names(df_list) <- c("E_rel")
      }
    }
    # Filtre per generar només els grafics de FARDEEP
    if (method=="FARDEEP"){
      # Filtre per fer només el gràfic de deconvolució relativa
      if (type=="rel"){
        F_rel <- FARDEEP_rel_NF(matrix, sig.matrix, results_dir, height, name)
        df_list <- list(F_rel)
        names(df_list) <- c("F_rel")
      }
      # Filtre per fer només el gràfic de deconvolució absoluta
      if (type=="abs"){
        F_abs <- FARDEEP_abs_NF(matrix, sig.matrix, results_dir, height, name)
        df_list <- list(F_abs)
        names(df_list) <- c("F_abs")
      }
      # Filtre per fer els grafics de deconvolució relativa i absoluta
      if (type=="abs_rel"){
        F_rel <- FARDEEP_rel_NF(matrix, sig.matrix, results_dir, height, name)
        F_abs <- FARDEEP_abs_NF(matrix, sig.matrix, results_dir, height, name)
        df_list <- list(F_rel, F_abs)
        names(df_list) <- c("F_rel", "F_abs")
      }
    }
  }
  # Primera selecció, entra SI es volen generar grafics fraccionats
  if (fraccionate.samples==TRUE){
    # Filtre per fer tots els metodes
    if (method=="ALL"){
      # Filtre per fer només els gràfics de deconvolució relativa
      if (type=="rel"){
        C_rel <- CIBERSORT_rel_F(matrix, sig.matrix, fractions, results_dir, height, name)
      }
      # Filtre per fer només els gràfics de deconvolució absoluta
      if (type=="abs"){
        C_abs <- CIBERSORT_abs_F(matrix, sig.matrix, fractions, results_dir, height, name)
      }
      # Filtre per fer només els gràfics de deconvolució relativa i absoluta
      if (type=="abs_rel"){
        C_rel <- CIBERSORT_rel_F(matrix, sig.matrix, fractions, results_dir, height, name)
        C_abs <- CIBERSORT_abs_F(matrix, sig.matrix, fractions, results_dir, height, name)
      }
      # Filtre per fer només els gràfics de deconvolució relativa
      if (type=="rel"){
        E_rel <- EPIC_rel_F(matrix, sig.matrix, fractions, results_dir, height, name)
      }
      # Filtre per fer només els gràfics de deconvolució absoluta, en aquest cas un warning perque no hi ha EPIC absolut
      if (type=="abs"){
        warning("Warning: EPIC only has relative type.")
      }
      # Filtre per fer només els gràfics de deconvolució relativa i absoluta, en aquest cas només els relatius
      if (type=="abs_rel"){
        E_rel <- EPIC_rel_F(matrix, sig.matrix, fractions, results_dir, height, name)
        warning("Warning: EPIC only has relative type.")
      }
      # Filtre per fer només els gràfics de deconvolució relativa
      if (type=="rel"){
        F_rel <- FARDEEP_rel_F(matrix, sig.matrix, fractions, results_dir, height, name)
        df_list <- list(C_rel, E_rel, F_rel)
        names(df_list) <- c("C_rel", "E_rel", "F_rel")
      }
      # Filtre per fer només els gràfics de deconvolució absoluta
      if (type=="abs"){
        F_abs <- FARDEEP_abs_F(matrix, sig.matrix, fractions, results_dir, height, name)
        df_list <- list(C_abs, F_abs)
        names(df_list) <- c("C_abs", "F_abs")
      }
      # Filtre per fer només els gràfics de deconvolució relativa i absoluta
      if (type=="abs_rel"){
        F_rel <- FARDEEP_rel_F(matrix, sig.matrix, fractions, results_dir, height, name)
        F_abs <- FARDEEP_abs_F(matrix, sig.matrix, fractions, results_dir, height, name)
        df_list <-list(C_rel, C_abs, E_rel, F_rel, F_abs)
        names(df_list) <- c("C_rel", "C_abs", "E_rel", "F_rel", "F_abs")
      }
    }
    # Filtre per fer els grafics amb el metode CIBERSORT
    if (method=="CIBERSORT"){
      # Filtre per fer només els gràfics de deconvolució relativa
      if (type=="rel"){
        C_rel <- CIBERSORT_rel_F(matrix, sig.matrix, fractions, results_dir, height, name)
        df_list <- list(C_rel)
        names(df_list) <- c("C_rel")
      }
      # Filtre per fer només els gràfics de deconvolució absoluta
      if (type=="abs"){
        C_abs <- CIBERSORT_abs_F(matrix, sig.matrix, fractions, results_dir, height, name)
        df_list <- list(C_abs)
        names(df_list) <- c("C_abs")
      }
      # Filtre per fer només els gràfics de deconvolució relativa i absoluta
      if (type=="abs_rel"){
        C_rel <- CIBERSORT_rel_F(matrix, sig.matrix, fractions, results_dir, height, name)
        C_abs <- CIBERSORT_abs_F(matrix, sig.matrix, fractions, results_dir, height, name)
        df_list <- list(C_rel, C_abs)
        names(df_list) <- c("C_rel", "C_abs")
      }
    }
    # Filtre per fer els grafics amb el metode EPIC
    if (method=="EPIC"){
      # Filtre per fer només els gràfics de deconvolució relativa
      if (type=="rel"){
        E_rel <- EPIC_rel_F(matrix, sig.matrix, fractions, results_dir, height, name)
        df_list <- list(E_rel)
        names(df_list) <- c("E_rel")
      }
      # Filtre per fer només els gràfics de deconvolució absoluta, en aquest cas un warning
      if (type=="abs"){
        warning("Warning: EPIC only has relative type.")
        df_list <- list()
      }
      # Filtre per fer només els gràfics de deconvolució relativa i el warning per absoluta
      if (type=="abs_rel"){
        E_rel <- EPIC_rel_F(matrix, sig.matrix, fractions, results_dir, height, name)
        warning("Warning: EPIC only has relative type.")
        df_list <- list(E_rel)
        names(df_list) <- c("E_rel")
      }
    }
    # Filtre per fer els grafics amb el metode FARDEEP
    if (method=="FARDEEP"){
      # Filtre per fer només els gràfics de deconvolució relativa
      if (type=="rel"){
        F_rel <- FARDEEP_rel_F(matrix, sig.matrix, fractions, results_dir, height, name)
        df_list <- list(F_rel)
        names(df_list) <- c("F_rel")
      }
      # Filtre per fer només els gràfics de deconvolució absoluta
      if (type=="abs"){
        F_abs <- FARDEEP_abs_F(matrix, sig.matrix, fractions, results_dir, height, name)
        df_list <- list(F_abs)
        names(df_list) <- c("F_abs")
      }
      # Filtre per fer només els gràfics de deconvolució relativa i absoluta
      if (type=="abs_rel"){
        F_rel <- FARDEEP_rel_F(matrix, sig.matrix, fractions, results_dir, height, name)
        F_abs <- FARDEEP_abs_F(matrix, sig.matrix, fractions, results_dir, height, name)
        df_list <- list(F_rel, F_abs)
        names(df_list) <- c("F_rel", "F_abs")
      }
    }
  }
  return(df_list)
}

