
##' @name Box_Deconv
##' @alias Box_Deconv
##' @title Box_Deconv
##' Box_Deconv
##'
##' Function that generates a boxplot comparing two conditions for each cell type from a df and a condition vector.
##' @param data Data Frame resulting from the deconvolution.
##' @param cond Vector that assigns a condition to each sample of the data frame.
##' @param results_dir Directory to which the graphics will be generated.
##' @param f_name Name you want to give to the file, by default it is "Deconvolution Boxplot".
##' @return The function returns the boxplots graph in the specified directory.
##' @author Nidia Barco Armengol
##' @export
##'
##' @examples
##' Box_Deconv(data= C_abs, cond= fractions, results_dir = "/bicoh/nidia/Deconv/Brain_MMartin/test", f_name= "Box_CIBERSORT_abs")

Box_Deconv <- function(data, cond, results_dir, f_name= "Deconvolution Boxplot"){

  require(ggplot2)
  require(ggpubr)
  require(gplots)
  require(data.table)
  require(ggsignif)
  require(tidyverse)

  # Manipulació del df de deconvolució, perque tingui l'estructura per generar els boxplots
  data$cell_type <- NULL
  data <- t(data)
  test_m <- melt(as.data.table(data), variable.factor = FALSE)
  condition <- c(rep(cond, times=length(colnames(data))))
  df_m <- cbind(test_m, condition)

  # Assignació del nom per al gràfic, per defecte "Deconvolution Boxplot"
  file <- c(f_name, ".png")
  file <- paste(file, collapse = "")

  # Si hi ha dues condicions es genera el grafic partit per cell type i amb el P_value per condició i apareix en cada caixa
  # Si hi ha mes de dues condicions es genera el grafic partit per cell type sense els P_value, els P_value per cada cell type surten com a missatges, i retorna un df amb info dels P_values i els nivells de significancia per condició
  if (length(unique(cond)) == 2){
    p <- ggplot(data = df_m, aes(x=variable, y=value)) +
      geom_boxplot(aes(fill=condition)) + stat_compare_means(aes(group=condition), method = "t.test")
    p + facet_wrap( ~ variable, scales="free")
    ggsave(filename=paste(results_dir,file, sep = "/"),dpi=300,width =10, height = 10)
  }else{

    list_P.value <- list()
    unics <- unique(df_m$variable)
    for (i in unics){
      positions <- c()
      for (x in 1:length(df_m$variable)){
        element <- df_m$variable[x]
        if (element==i){
          positions <- c(positions, x)
          list_P.value[[i]] <- positions
        }
      }
    }
    for (x in 1:length(list_P.value)){
      l <- list_P.value[[x]]
      v <- c()
      for (e in l){
        v <- c(v, df_m$value[e])
      }
      list_P.value[[x]] <- v
    }
    for (v in 1:length(list_P.value)){
      l <- list_P.value[[v]]
      p <- t.test(l, mu=0)$p.value
      list_P.value[[v]] <- p
    }
    list_P.value <- data_frame(list_P.value)
    list_P.value <- cbind(unics, list_P.value)
    for (e in 1:length(rownames(list_P.value))){
      f <- list_P.value[e,1]
      f2 <- list_P.value[e,2]
      message("P.Value -> ",unlist(f),"=",unlist(f2))
    }

    p <- ggplot(data = df_m, aes(x=variable, y=value)) +
      geom_boxplot(aes(fill=condition))
    p + facet_wrap( ~ variable, scales="free")
    ggsave(filename=paste(results_dir,file, sep = "/"),dpi=300,width =10, height = 10)
    anno_df = compare_means(value ~ condition, group.by = "variable", data = df_m, method = "t.test")
    return(anno_df)
  }
}
