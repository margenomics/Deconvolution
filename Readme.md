Readme
================

# Introduction

The Deconvolution package contains different functions, specific to
generate absolute or relative deconvolution plots, with different
methods, which also return data frames with the result of the
deconvolution, or a list with data frames in the case of the deconv
function, these df can then be used with another function of the
package, the Box\_Deconv function to generate a boxplot that compares
two conditions in the different cell types, and does so from the df
resulting from the deconvolution.

# Functions

| Function           | Description                                                                                                                                                         | Return                                    |
|--------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------|
| CIBERSORT\_abs\_F  | Function that generates the absolute deconvolution plot with the CIBERSORT method, with the fractional plot, and returns a df resulting from the deconvolution.     | Deconvolution Graphs and Deconvolution df |
| CIBERSORT\_abs\_NF | Function that generates the absolute deconvolution plot with the CIBERSORT method, with the non-fractional plot, and returns a df resulting from the deconvolution. | Deconvolution Graph and df                |
| CIBERSORT\_rel\_F  | Function that generates the relative deconvolution plot with the CIBERSORT method, with the fractional plot, and returns a df resulting from the deconvolution.     | Deconvolution Graphs and Deconvolution df |
| CIBERSORT\_rel\_NF | Function that generates the relative deconvolution plot with the CIBERSORT method, with the non-fractional plot, and returns a df resulting from the deconvolution. | Deconvolution Graph and df                |
| EPIC\_rel\_F       | Function that generates the relative deconvolution plot with the EPIC method, with the fractional plot, and returns a df resulting from the deconvolution.          | Deconvolution Graphs and Deconvolution df |
| EPIC\_rel\_NF      | Function that generates the relative deconvolution plot with the EPIC method, with the non-fractional plot, and returns a df resulting from the deconvolution.      | Deconvolution Graph and df                |
| FARDEEP\_abs\_F    | Function that generates the absolute deconvolution plot with the FARDEEP method, with the fractional plot, and returns a df resulting from the deconvolution.       | Deconvolution Graphs and Deconvolution df |
| FARDEEP\_abs\_NF   | Function that generates the absolute deconvolution plot with the FARDEEP method, with the non-fractional plot, and returns a df resulting from the deconvolution.   | Deconvolution Graph and df                |
| FARDEEP\_rel\_F    | Function that generates the relative deconvolution plot with the FARDEEP method, with the fractional plot, and returns a df resulting from the deconvolution.       | Deconvolution Graphs and Deconvolution df |
| FARDEEP\_rel\_NF   | Function that generates the relative deconvolution plot with the FARDEEP method, with the non-fractional plot, and returns a df resulting from the deconvolution.   | Deconvolution Graph and df                |
| Deconv             | Function that generates massive deconvolution plots with the CIBERSORT, EPIC and FARDEEP methods, and returns a list of df’s resulting from the deconvolutions.     | Deconvolution Graphs and df’s list        |
| Box\_Deconv        | Function that generates a boxplot comparing two conditions for each cell type from a df and a condition vector.                                                     | Deconvolution Boxplot                     |

# Installation

``` r
library(devtools)
devtools::install_github("margenomics/Deconvolution")
library(Deconvolution)
```

# Usage

``` r
# CIBERSORT_abs_F
  C_abs <- CIBERSORT_abs_F("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", "Monaco-ABIS_Immune11_Microarray.txt", f, results_dir = "/bicoh/nidia/Deconv/Brain_MMartin/test", height= 8, name, "MONACO")

# CIBERSORT_abs_NF
  C_abs <- CIBERSORT_abs_NF("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", "Monaco-ABIS_Immune11_Microarray.txt", results_dir = "/bicoh/nidia/Deconv/Brain_MMartin/test", height= 8, name= "MONACO")

# CIBERSORT_rel_F
  C_rel <- CIBERSORT_rel_F("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", "Monaco-ABIS_Immune11_Microarray.txt", f, results_dir = "/bicoh/nidia/Deconv/Brain_MMartin/test", height= 8, name= "MONACO")

# CIBERSORT_rel_NF
  C_rel <- CIBERSORT_rel_NF("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", "Monaco-ABIS_Immune11_Microarray.txt", results_dir = "/bicoh/nidia/Deconv/Brain_MMartin/test", height= 8, name= "MONACO")

# EPIC_rel_F
  E_rel <- EPIC_rel_F("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", "Monaco-ABIS_Immune11_Microarray.txt", f, results_dir = "/bicoh/nidia/Deconv/Brain_MMartin/test", height= 8, name= "MONACO")

# EPIC_rel_NF
  E_rel <- EPIC_rel_NF("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", "Monaco-ABIS_Immune11_Microarray.txt", results_dir = "/bicoh/nidia/Deconv/Brain_MMartin/test", height= 8, name= "MONACO")

# FARDEEP_abs_F
  F_abs <- FARDEEP_abs_F("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", "Monaco-ABIS_Immune11_Microarray.txt", f, results_dir = "/bicoh/nidia/Deconv/Brain_MMartin/test", height= 8, name= "MONACO")

# FARDEEP_abs_NF
  F_abs <- FARDEEP_abs_NF("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", "Monaco-ABIS_Immune11_Microarray.txt", results_dir = "/bicoh/nidia/Deconv/Brain_MMartin/test", height= 8, name= "MONACO")

# FARDEEP_rel_F
  F_rel <- FARDEEP_rel_F("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", "Monaco-ABIS_Immune11_Microarray.txt", f, results_dir = "/bicoh/nidia/Deconv/Brain_MMartin/test", height= 8, name= "MONACO")

# FARDEEP_rel_NF
  F_rel <- FARDEEP_rel_NF("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", "Monaco-ABIS_Immune11_Microarray.txt", results_dir = "/bicoh/nidia/Deconv/Brain_MMartin/test", height= 8, name= "MONACO")

# Deconv
  L_rel <- deconv("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", "Monaco-ABIS_Immune11_Microarray.txt", "ALL", "rel", fraccionate.samples=TRUE, fractions= f, results_dir = "/bicoh/nidia/Deconv/Brain_MMartin/test", height= 8, name= "MONACO")

# Box_deconv
  P.value_df <- Box_Deconv(data= C_abs, cond= fractions, results_dir = "/bicoh/nidia/Deconv/Brain_MMartin/test", f_name= "Box_CIBERSORT_abs")
```

# Examples

``` r
matrix <- "Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt"
sig.matrix <- "Monaco-ABIS_Immune11_Microarray.txt"
fractions <- c(rep("A", each=18), rep("B", each=43), rep("C", each=36))
results_dir <- "/bicoh/nidia/Deconv/Brain_MMartin/test"
```

## CIBERSORT\_abs\_NF

### Itroduce:

``` r
C_abs <- CIBERSORT_abs_NF("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", "Monaco-ABIS_Immune11_Microarray.txt", results_dir = "/bicoh/nidia/Deconv/Brain_MMartin/test")
```

### Return:

![Image
text](/bicoh/nidia/Deconvolution/Images/CIBERSORT.abs_Monaco-ABIS_Immune11_Microarray_plot.png)

``` r
# Show the first 6 columns of the deconvolution df.
knitr::kable(C_abs[,1:6 ])
```

|                | Microgl.f.1\_PAG | Microgl.f.P1\_SC | Microgl.f.P1\_SS | Microgl.f.13\_PAG | Microgl.f.P13\_SC | Microgl.f.P13\_SS |
|:---------------|-----------------:|-----------------:|-----------------:|------------------:|------------------:|------------------:|
| T Naive        |        0.0361892 |        0.0000000 |        0.0039591 |         0.0168165 |         0.0131280 |         0.0208645 |
| T Memory       |        0.0000000 |        0.0067183 |        0.0000000 |         0.0000000 |         0.0011728 |         0.0000000 |
| B Naive        |        0.0000000 |        0.0053894 |        0.0000000 |         0.0000000 |         0.0441601 |         0.0000000 |
| B Memory       |        0.0000000 |        0.0000000 |        0.0040544 |         0.0000000 |         0.0000000 |         0.0000000 |
| Plasmablasts   |        0.0041236 |        0.0037040 |        0.0044582 |         0.0009182 |         0.0046624 |         0.0051573 |
| NK             |        0.0000000 |        0.0000000 |        0.0000000 |         0.0035398 |         0.0000000 |         0.0007784 |
| pDCs           |        0.1506466 |        0.1657560 |        0.1218028 |         0.1211761 |         0.1983397 |         0.1274756 |
| Neutrophils LD |        0.0000000 |        0.0000000 |        0.0000000 |         0.0000000 |         0.0000000 |         0.0000000 |
| Basophils LD   |        0.0101286 |        0.0148763 |        0.0004157 |         0.0000000 |         0.0000000 |         0.0017440 |
| mDCs           |        0.0115032 |        0.0220572 |        0.0094786 |         0.0165953 |         0.0299044 |         0.0089637 |
| Monocytes      |        0.1231687 |        0.0778489 |        0.1902689 |         0.0984570 |         0.0870736 |         0.1661126 |

## CIBERSORT\_rel\_F

### Itroduce:

``` r
C_rel <- CIBERSORT_rel_F("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", "Monaco-ABIS_Immune11_Microarray.txt", f, results_dir = "/bicoh/nidia/Deconv/Brain_MMartin/test", height= 8, name= "MONACO")
```

### Return:

![Image text](/bicoh/nidia/Deconvolution/Images/C_rel_fraccions.png)

``` r
# Show the first 6 columns of the deconvolution df.
knitr::kable(C_res[,1:6 ])
```

|                | Microgl.f.1\_PAG | Microgl.f.P1\_SC | Microgl.f.P1\_SS | Microgl.f.13\_PAG | Microgl.f.P13\_SC | Microgl.f.P13\_SS |
|:---------------|-----------------:|-----------------:|-----------------:|------------------:|------------------:|------------------:|
| T Naive        |        0.1077831 |        0.0000000 |        0.0118381 |         0.0653062 |         0.0346896 |         0.0630164 |
| T Memory       |        0.0000000 |        0.0226702 |        0.0000000 |         0.0000000 |         0.0030991 |         0.0000000 |
| B Naive        |        0.0000000 |        0.0181858 |        0.0000000 |         0.0000000 |         0.1166894 |         0.0000000 |
| B Memory       |        0.0000000 |        0.0000000 |        0.0121229 |         0.0000000 |         0.0000000 |         0.0000000 |
| Plasmablasts   |        0.0122813 |        0.0124986 |        0.0133304 |         0.0035658 |         0.0123201 |         0.0155764 |
| NK             |        0.0000000 |        0.0000000 |        0.0000000 |         0.0137468 |         0.0000000 |         0.0023508 |
| pDCs           |        0.4486736 |        0.5593251 |        0.3642019 |         0.4705813 |         0.5240969 |         0.3850110 |
| Neutrophils LD |        0.0000000 |        0.0000000 |        0.0000000 |         0.0000000 |         0.0000000 |         0.0000000 |
| Basophils LD   |        0.0301663 |        0.0501983 |        0.0012429 |         0.0000000 |         0.0000000 |         0.0052673 |
| mDCs           |        0.0342601 |        0.0744296 |        0.0283418 |         0.0644472 |         0.0790199 |         0.0270728 |
