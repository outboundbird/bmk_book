#' ---
#' title:
#' subtitle: 'SAR: sar , Study: study'
#' author:  Siying Huang (E0482362), Biomarker statistics team
#' date: 'created: 2023-07-28 , updated (`r Sys.Date()`)'
#' always_allow_html: true
#' output:
#'   html_document:
#'     css: theme.css
#'     code_folding: "hide"
#'     toc: yes
#'     toc_float:
#'       collapse: no
#' ---
#+ setup, include = FALSE
knitr::opts_knit$set(root_dir = "/mnt/c/Users/e0482362/Work/bookdown/bookdown/src")
knitr::opts_chunk$set(echo = T, comment = "", message = F, warning = F, error = F)
options(width = 100)
#' # Quantitative PCR (qPCR)
#' qPCR is another method to quantify the gene expressions. Unlike other high-throughput methods, qPCR usually is not genomewide covered and targets only a few genes at a time.
#'
#' Here is a youtube tutorial on [Overview of qPCR](https://youtu.be/1kvy17ugI4w).
#' [tutorial on data analysis of qPCR](https://youtu.be/GQOnX1-SUrI)
#'
#' ## Quantification of gene expression levels
#' ![Fig. courtesy of BIO-RAD tutorial](https://www.bio-rad.com/webroot/web/images/lsr/solutions/technologies/gene_expression/qPCR_real-time_PCR/technology_detail/real-time-pcr-amplification-plot-qpcr-real-time-pcr.jpg)
#'
#' In general there are two basic typyes of RT-qPCR quantification analysis:
#'
#' - Absolute quantification. To answer the quation of 'How many?', such as chromosome or gene copy number determindation and viral load measurements
#' - Relative quantification. to compare the changes in gene expression, to answer the question of  what is the fold difference? Usually involves multiple genes.
#' Here I elaborate more on the relative quantification method.
#' ## Delta-delta method $2^{-\Delta \Delta t}$

#+ libs
library(here)
library(learnr)

#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/061_qPCR.R', output_dir = 'output')
# knitr::spin('src/061_qPCR.R', format = 'Rmd', knit = FALSE)