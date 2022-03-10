Readme
================
Nidia Barco Armengol

# Introduction

The Deconvolution package contains different functions, specific to
generate absolute or relative deconvolution plots, with different
methods, which also return data frames with the result of the
deconvolution, or a list with data frames in the case of the deconv
function, these df can then be used with another function of the
package, the Box\_Deconv function to generate a boxplot that compares
two conditions in the different cell types, and does so from the df
resulting from the deconvolution.

# Methods

| Deconvolution method                                           | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
|----------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| CIBERSORT                                                      | Versatile computational method for quantifying cell fractions from bulk tissue gene expression profiles (GEPs). Combining support vector regression with prior knowledge of expression profiles from purified leukocyte subsets,this methode can accurately estimate the immune composition of a tumor biopsy. <https://doi.org/10.1007/978-1-4939-7493-1_12>                                                                                                                                      |
| EPIC (Estimating the Proportions of Immune and Cancer cells)   | includes RNA-seq-based gene expression reference profiles from immune cells and other nonmalignant cell types found in tumors, can additionally manage user-defined gene expression reference profiles, include the ability to account for an uncharacterized cell type. <https://doi.org/10.1007/978-1-0716-0327-7_17>                                                                                                                                                                            |
| FARDEEP (Fast And Robust DEconvolution of Expression Profiles) | This methode enumerate immune cell subsets from whole tumor tissue samples. To reduce noise in the tumor gene expression datasets, utilizes an adaptive least trimmed square to automatically detect and remove outliers before estimating the cell compositions. We show that FARDEEP is less susceptible to outliers and returns a better estimation of coefficients than the existing methods with both numerical simulations and real datasets. <https://doi.org/10.1371/journal.pcbi.1006976> |

# Signature matrix

| Signature matrix                  | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
|-----------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| LM22                              | Is a signature matrix file consisting of 547 genes that accurately distinguish 22 mature human hematopoietic populations isolated from peripheral blood or in vitro culture conditions, including seven T cell types, naïve and memory B cells, plasma cells, NK cells, and myeloid subsets. <https://doi.org/10.1007/978-1-4939-7493-1_12>                                                                                                                                                                                                                 |
| Monaco-ABIS\_Immune11\_Microarray | This matrix was constructed from a set of predictor variables (cell types), normalized by mRNA abundance with scaling factors, and filtered by extracting genes with very high or very low expression and also genes not expressed with the STAR method, including eleven cell tipes, T.Naive, T.Memory, B.Naive, B.Memory, Plasmablasts, NK, pDCs, Neutrophils.LD, Basophils.LD, mDCs, Monocytes. <https://doi.org/10.1016/j.celrep.2019.01.041>                                                                                                           |
| Wang2020\_signature.matrix        | To construct the signature matrix of brain data, adult scRNA-seq data from Darmanis et al. (2015) were analyzed by grouping characterized human brain cells into seven cell types \[astrocytes, endothelia, microglia, oligodendrocytes (Oligo), oligodendrocyte precursor cells (OPCs), inhibitory and excitatory neurons\] and selecting the top 50 marker genes for each cell type using SC3 (Kiselev et al., 2017). Finally, the expression of each marker gene in cells of the same type was averaged. <https://doi.org/10.1093/bioinformatics/btz619> |

# Functions

| Function           | Description                                                                                                                                                                                             | Return                                    |
|--------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------|
| CIBERSORT\_abs\_F  | Function that generates the absolute deconvolution plot with the CIBERSORT method, with the fractional plot, and returns a df resulting from the deconvolution.                                         | Deconvolution Graphs and Deconvolution df |
| CIBERSORT\_abs\_NF | Function that generates the absolute deconvolution plot with the CIBERSORT method, with the non-fractional plot, and returns a df resulting from the deconvolution.                                     | Deconvolution Graph and df                |
| CIBERSORT\_rel\_F  | Function that generates the relative deconvolution plot with the CIBERSORT method, with the fractional plot, and returns a df resulting from the deconvolution.                                         | Deconvolution Graphs and Deconvolution df |
| CIBERSORT\_rel\_NF | Function that generates the relative deconvolution plot with the CIBERSORT method, with the non-fractional plot, and returns a df resulting from the deconvolution.                                     | Deconvolution Graph and df                |
| EPIC\_rel\_F       | Function that generates the relative deconvolution plot with the EPIC method, with the fractional plot, and returns a df resulting from the deconvolution.                                              | Deconvolution Graphs and Deconvolution df |
| EPIC\_rel\_NF      | Function that generates the relative deconvolution plot with the EPIC method, with the non-fractional plot, and returns a df resulting from the deconvolution.                                          | Deconvolution Graph and df                |
| FARDEEP\_abs\_F    | Function that generates the absolute deconvolution plot with the FARDEEP method, with the fractional plot, and returns a df resulting from the deconvolution.                                           | Deconvolution Graphs and Deconvolution df |
| FARDEEP\_abs\_NF   | Function that generates the absolute deconvolution plot with the FARDEEP method, with the non-fractional plot, and returns a df resulting from the deconvolution.                                       | Deconvolution Graph and df                |
| FARDEEP\_rel\_F    | Function that generates the relative deconvolution plot with the FARDEEP method, with the fractional plot, and returns a df resulting from the deconvolution.                                           | Deconvolution Graphs and Deconvolution df |
| FARDEEP\_rel\_NF   | Function that generates the relative deconvolution plot with the FARDEEP method, with the non-fractional plot, and returns a df resulting from the deconvolution.                                       | Deconvolution Graph and df                |
| Deconv             | Function that generates massive deconvolution plots with the CIBERSORT, EPIC and FARDEEP methods, and returns a list of df’s resulting from the deconvolutions.                                         | Deconvolution Graphs and df’s list        |
| Box\_Deconv        | Function that generates a boxplot comparing conditions for each cell type from a df and a condition vector, and in the case of more than two conditions, it also returns a df with P.Value information. | Deconvolution Boxplot                     |

# Installation

``` r
library(devtools)
devtools::install_github("margenomics/Deconvolution")
library(Deconvolution)
```

# Usage

``` r
# CIBERSORT_abs_F
  C_abs <- CIBERSORT_abs_F("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", 
  "Wang2020_signature.matrix.txt", f, results_dir = "/bicoh/nidia/Deconv/Brain/test", height= 8, name, "MONACO")

# CIBERSORT_abs_NF
  C_abs <- CIBERSORT_abs_NF("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", 
  "Wang2020_signature.matrix.txt", results_dir = "/bicoh/nidia/Deconv/Brain/test", height= 8, name= "MONACO")

# CIBERSORT_rel_F
  C_rel <- CIBERSORT_rel_F("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", 
  "Wang2020_signature.matrix.txt", f, results_dir = "/bicoh/nidia/Deconv/Brain/test", height= 8, name= "MONACO")

# CIBERSORT_rel_NF
  C_rel <- CIBERSORT_rel_NF("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", 
  "Wang2020_signature.matrix.txt", results_dir = "/bicoh/nidia/Deconv/Brain/test", height= 8, name= "MONACO")

# EPIC_rel_F
  E_rel <- EPIC_rel_F("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", 
  "Wang2020_signature.matrix.txt", f, results_dir = "/bicoh/nidia/Deconv/Brain/test", height= 8, name= "MONACO")

# EPIC_rel_NF
  E_rel <- EPIC_rel_NF("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", 
  "Wang2020_signature.matrix.txt", results_dir = "/bicoh/nidia/Deconv/Brain/test", height= 8, name= "MONACO")

# FARDEEP_abs_F
  F_abs <- FARDEEP_abs_F("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", 
  "Wang2020_signature.matrix.txt", f, results_dir = "/bicoh/nidia/Deconv/Brain/test", height= 8, name= "MONACO")

# FARDEEP_abs_NF
  F_abs <- FARDEEP_abs_NF("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", 
  "Wang2020_signature.matrix.txt", results_dir = "/bicoh/nidia/Deconv/Brain/test", height= 8, name= "MONACO")

# FARDEEP_rel_F
  F_rel <- FARDEEP_rel_F("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", 
  "Wang2020_signature.matrix.txt", f, results_dir = "/bicoh/nidia/Deconv/Brain/test", height= 8, name= "MONACO")

# FARDEEP_rel_NF
  F_rel <- FARDEEP_rel_NF("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", 
  "Wang2020_signature.matrix.txt", results_dir = "/bicoh/nidia/Deconv/Brain/test", height= 8, name= "MONACO")

# Deconv
  L_rel <- deconv("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", 
  "Wang2020_signature.matrix.txt", "ALL", "rel", fraccionate.samples=TRUE, 
  fractions= f, results_dir = "/bicoh/nidia/Deconv/Brain/test", height= 8, name= "MONACO")

# Box_deconv
  P.value_df <- Box_Deconv(data= C_abs, cond= fractions, 
  results_dir = "/bicoh/nidia/Deconv/Brain/test", f_name= "Box_CIBERSORT_abs")
```

# Examples

``` r
matrix <- "Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt"
sig.matrix <- "Wang2020_signature.matrix.txt"
fractions <- c(rep("A", each=18), rep("B", each=43), rep("C", each=36))
results_dir <- "/bicoh/nidia/Deconv/Brain/test"
```

## CIBERSORT\_abs\_NF

### Itroduce:

``` r
# The CIBERSORT_abs_NF function generates the absolute deconvolution graph with the cibersort method, and returns the df resulting from the deconvolution. 

C_abs <- CIBERSORT_abs_NF("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", "Wang2020_signature.matrix.txt", results_dir = "/bicoh/nidia/Deconv/Brain/test")
```

### Return:

![Image
text](/Images/CIBERSORT.abs_Wang2020_signature.matrix__plot.png)

``` r
# Show the first 6 columns of the deconvolution df.
knitr::kable(C_abs[,1:6 ])
```

|             | Microgl.f.1\_PAG | Microgl.f.P1\_SC | Microgl.f.P1\_SS | Microgl.f.13\_PAG | Microgl.f.P13\_SC | Microgl.f.P13\_SS |
|:------------|-----------------:|-----------------:|-----------------:|------------------:|------------------:|------------------:|
| Astrocyte   |        0.0095227 |        0.0090841 |        0.0041076 |         0.0177518 |         0.0100331 |         0.0034658 |
| Endothelial |        0.0049448 |        0.0034026 |        0.0019898 |         0.0022893 |         0.0081100 |         0.0010536 |
| Microglia   |        0.0000000 |        0.0026684 |        0.0000000 |         0.0000000 |         0.0006046 |         0.0000000 |
| Excitatory  |        0.0068483 |        0.0050919 |        0.0145409 |         0.0103304 |         0.0099650 |         0.0187261 |
| Inhibitory  |        0.0125782 |        0.0040904 |        0.0000000 |         0.0192275 |         0.0009906 |         0.0000000 |
| Oligo       |        0.0071585 |        0.0127354 |        0.0020316 |         0.0091130 |         0.0124069 |         0.0000000 |

## CIBERSORT\_rel\_F

### Itroduce:

``` r
# The CIBERSORT_rel_F function generates the relative deconvolution graphs, as many as fractions in the fractions vector with the cibersort method, and returns the resulting deconvolution df. 

C_rel <- CIBERSORT_rel_F("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", "Wang2020_signature.matrix.txt", f, results_dir = "/bicoh/nidia/Deconv/Brain/test", height= 8, name= "MONACO")
```

### Return:

![Image text](/Images/C_rel_fractions.png)

``` r
# Show the first 6 columns of the deconvolution df.
knitr::kable(C_rel[,1:6 ])
```

|             | Microgl.f.1\_PAG | Microgl.f.P1\_SC | Microgl.f.P1\_SS | Microgl.f.13\_PAG | Microgl.f.P13\_SC | Microgl.f.P13\_SS |
|:------------|-----------------:|-----------------:|-----------------:|------------------:|------------------:|------------------:|
| Astrocyte   |        0.2319633 |        0.2450339 |        0.1811911 |         0.3023533 |         0.2382573 |         0.1490946 |
| Endothelial |        0.1204505 |        0.0917826 |        0.0877717 |         0.0389926 |         0.1925898 |         0.0453257 |
| Microglia   |        0.0000000 |        0.0719781 |        0.0000000 |         0.0000000 |         0.0143573 |         0.0000000 |
| Excitatory  |        0.1668184 |        0.1373495 |        0.6414209 |         0.1759511 |         0.2366410 |         0.8055797 |
| Inhibitory  |        0.3063941 |        0.1103333 |        0.0000000 |         0.3274878 |         0.0235249 |         0.0000000 |

## Deconv

### Introduce

``` r
# The deconv function can generate any deconvolution graph, fractional or not and with any of the three methods contained in the package (CIBERSORT, EPIC or FARDEEP), and returns a list of the df's resulting from the deconvolution with each method. 

L_rel <- deconv("Counts.HUMAN_nodup_Sal_Veh_IPs_INPUTs.txt", "Wang2020_signature.matrix.txt", "ALL", "rel", fraccionate.samples=FALSE, results_dir = "/bicoh/nidia/Deconv/Brain/test")
```

### Return

![Image text](/Images/Ex_deconv.png)

``` r
# l_df <- list(C_rel, E_rel, F_rel)
knitr::kable(l_df[[1]][,1:6])
```

|             | Microgl.f.1\_PAG | Microgl.f.P1\_SC | Microgl.f.P1\_SS | Microgl.f.13\_PAG | Microgl.f.P13\_SC | Microgl.f.P13\_SS |
|:------------|-----------------:|-----------------:|-----------------:|------------------:|------------------:|------------------:|
| Astrocyte   |        0.2319633 |        0.2450339 |        0.1811911 |         0.3023533 |         0.2382573 |         0.1490946 |
| Endothelial |        0.1204505 |        0.0917826 |        0.0877717 |         0.0389926 |         0.1925898 |         0.0453257 |
| Microglia   |        0.0000000 |        0.0719781 |        0.0000000 |         0.0000000 |         0.0143573 |         0.0000000 |
| Excitatory  |        0.1668184 |        0.1373495 |        0.6414209 |         0.1759511 |         0.2366410 |         0.8055797 |
| Inhibitory  |        0.3063941 |        0.1103333 |        0.0000000 |         0.3274878 |         0.0235249 |         0.0000000 |

``` r
knitr::kable(l_df[[2]][,1:6])
```

|             | Microgl.f.1\_PAG | Microgl.f.P1\_SC | Microgl.f.P1\_SS | Microgl.f.13\_PAG | Microgl.f.P13\_SC | Microgl.f.P13\_SS |
|:------------|-----------------:|-----------------:|-----------------:|------------------:|------------------:|------------------:|
| Astrocyte   |        0.1531795 |        0.1556428 |        0.2864268 |         0.2318391 |         0.2240954 |         0.1569012 |
| Endothelial |        0.4571354 |        0.1533534 |        0.4320599 |         0.4302060 |         0.2527313 |         0.5786744 |
| Microglia   |        0.1682455 |        0.1790524 |        0.2390194 |         0.1229658 |         0.1476324 |         0.2542964 |
| Excitatory  |        0.0000460 |        0.0000004 |        0.0000008 |         0.0000013 |         0.0036760 |         0.0000029 |
| Inhibitory  |        0.0000002 |        0.0000051 |        0.0165302 |         0.0044590 |         0.0765474 |         0.0054710 |
| Oligo       |        0.2213921 |        0.5119449 |        0.0259628 |         0.2105285 |         0.2953168 |         0.0046540 |
| otherCells  |        0.0000013 |        0.0000009 |        0.0000001 |         0.0000004 |         0.0000007 |         0.0000000 |

``` r
knitr::kable(l_df[[3]][,1:6])
```

|             | Microgl.f.1\_PAG | Microgl.f.P1\_SC | Microgl.f.P1\_SS | Microgl.f.13\_PAG | Microgl.f.P13\_SC | Microgl.f.P13\_SS |
|:------------|-----------------:|-----------------:|-----------------:|------------------:|------------------:|------------------:|
| Astrocyte   |        0.3470971 |        0.4039560 |        0.2704517 |         0.3044236 |         0.3303557 |         0.2351979 |
| Endothelial |        0.0000000 |        0.0000000 |        0.0000000 |         0.0000000 |         0.0000000 |         0.0000000 |
| Microglia   |        0.0817858 |        0.1573647 |        0.0607674 |         0.0274723 |         0.3272161 |         0.0821891 |
| Excitatory  |        0.2555004 |        0.2237835 |        0.5545924 |         0.3196995 |         0.1666687 |         0.6199324 |
| Inhibitory  |        0.2350374 |        0.0774297 |        0.0000000 |         0.3088119 |         0.0351121 |         0.0000000 |
| Oligo       |        0.0805793 |        0.1374661 |        0.1141886 |         0.0395927 |         0.1406473 |         0.0626806 |

## Box\_Deconv

### Introduce

``` r
# The Box_Deconv function generates a graph comparing for each cell type the conditions present in the list of samples. 

P.value_df <- Box_Deconv(data= C_abs, cond= fractions, results_dir = "/bicoh/nidia/Deconv/Brain/test")
```

### Return

#### With 2 conditions

If the samples have two conditions, the function returns a graph with
the P.Value of the conditions for each cell type in the graph.

![Image text](/Images/2_cond_ex.png)

#### With more than 2 conditions

If the samples have more than two conditions, the function returns a
graph comparing the conditions for each cell type, but the overall
P.Value for each cell type is displayed as a message, since the function
returns a df containing information of the P.value’s and the
significance level of the conditions per cell type.

![Image text](/Images/3_cond_ex.png)

``` r
knitr::kable(anno_df)
```

| variable    | .y.   | group1 | group2 |         p |   p.adj | p.format   | p.signif | method |
|:------------|:------|:-------|:-------|----------:|--------:|:-----------|:---------|:-------|
| Astrocyte   | value | A      | B      | 0.0000000 | 0.0e+00 | 2.5e-12    | \*\*\*\* | T-test |
| Astrocyte   | value | A      | C      | 0.0010372 | 7.3e-03 | 0.0010     | \*\*     | T-test |
| Astrocyte   | value | B      | C      | 0.0000000 | 0.0e+00 | &lt; 2e-16 | \*\*\*\* | T-test |
| Endothelial | value | A      | B      | 0.1218155 | 2.4e-01 | 0.1218     | ns       | T-test |
| Endothelial | value | A      | C      | 0.0000079 | 7.9e-05 | 7.9e-06    | \*\*\*\* | T-test |
| Endothelial | value | B      | C      | 0.0000000 | 4.0e-07 | 3.3e-08    | \*\*\*\* | T-test |
| Microglia   | value | A      | B      | 0.0089715 | 4.5e-02 | 0.0090     | \*\*     | T-test |
| Microglia   | value | A      | C      | 0.0089715 | 4.5e-02 | 0.0090     | \*\*     | T-test |
| Excitatory  | value | A      | B      | 0.0662205 | 2.0e-01 | 0.0662     | ns       | T-test |
| Excitatory  | value | A      | C      | 0.0000000 | 4.0e-07 | 3.4e-08    | \*\*\*\* | T-test |
| Excitatory  | value | B      | C      | 0.0000041 | 4.5e-05 | 4.1e-06    | \*\*\*\* | T-test |
| Inhibitory  | value | A      | B      | 0.1927551 | 2.4e-01 | 0.1928     | ns       | T-test |
| Inhibitory  | value | A      | C      | 0.0000126 | 1.1e-04 | 1.3e-05    | \*\*\*\* | T-test |
| Inhibitory  | value | B      | C      | 0.0000436 | 3.5e-04 | 4.4e-05    | \*\*\*\* | T-test |
| Oligo       | value | A      | B      | 0.0013822 | 8.3e-03 | 0.0014     | \*\*     | T-test |
| Oligo       | value | A      | C      | 0.0000000 | 0.0e+00 | &lt; 2e-16 | \*\*\*\* | T-test |
| Oligo       | value | B      | C      | 0.0000000 | 0.0e+00 | 5.3e-12    | \*\*\*\* | T-test |

# References

Hao, Y., Yan, M., Heath, B. R., Lei, Y. L., & Xie, Y. (2019). Fast and
robust deconvolution of tumor infiltrating lymphocyte from expression
profiles using least trimmed squares. PLoS computational biology, 15(5),
e1006976. <https://doi.org/10.1371/journal.pcbi.1006976>

Chen, B., Khodadoust, M. S., Liu, C. L., Newman, A. M., & Alizadeh, A.
A. (2018). Profiling Tumor Infiltrating Immune Cells with CIBERSORT.
Methods in molecular biology (Clifton, N.J.), 1711, 243–259.
<https://doi.org/10.1007/978-1-4939-7493-1_12>

Racle, J., & Gfeller, D. (2020). EPIC: A Tool to Estimate the
Proportions of Different Cell Types from Bulk Gene Expression Data.
Methods in molecular biology (Clifton, N.J.), 2120, 233–248.
<https://doi.org/10.1007/978-1-0716-0327-7_17>

Avila Cobos, F., Alquicira-Hernandez, J., Powell, J.E. et
al. Benchmarking of cell type deconvolution pipelines for
transcriptomics data. Nat Commun 11, 5650 (2020).
<https://doi.org/10.1038/s41467-020-19015-1>

Gregor Sturm, Francesca Finotello, Florent Petitprez, Jitao David Zhang,
Jan Baumbach, Wolf H Fridman, Markus List, Tatsiana Aneichyk,
Comprehensive evaluation of transcriptome-based cell-type quantification
methods for immuno-oncology, Bioinformatics, Volume 35, Issue 14, July
2019, Pages i436–i445, <https://doi.org/10.1093/bioinformatics/btz363>

Francisco Avila Cobos, Jo Vandesompele, Pieter Mestdagh, Katleen De
Preter, Computational deconvolution of transcriptomics data from mixed
cell populations, Bioinformatics, Volume 34, Issue 11, 01 June 2018,
Pages 1969–1979, <https://doi.org/10.1093/bioinformatics/bty019>

Finotello, F., Trajanoski, Z. Quantifying tumor-infiltrating immune
cells from transcriptomics data. Cancer Immunol Immunother 67, 1031–1040
(2018). <https://doi.org/10.1007/s00262-018-2150-z>

Monaco, G., Lee, B., Xu, W., Mustafah, S., Hwang, Y. Y., Carré, C.,
Burdin, N., Visan, L., Ceccarelli, M., Poidinger, M., Zippelius, A.,
Pedro de Magalhães, J., & Larbi, A. (2019). RNA-Seq Signatures
Normalized by mRNA Abundance Allow Absolute Deconvolution of Human
Immune Cell Types. Cell reports, 26(6), 1627–1640.e7.
<https://doi.org/10.1016/j.celrep.2019.01.041>

Jiebiao Wang, Bernie Devlin, Kathryn Roeder, Using multiple measurements
of tissue to estimate subject- and cell-type-specific gene expression,
Bioinformatics, Volume 36, Issue 3, 1 February 2020, Pages 782–788,
<https://doi.org/10.1093/bioinformatics/btz619>
