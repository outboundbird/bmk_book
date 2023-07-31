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


#+ libs
library(here)
library(rstatix)
library(dplyr)
#+ io
qpcr <- readRDS(file.path(here(), "data/qPCR.rds"))
bsl <- qpcr %>%
  filter(Day == "D1")
head(bsl)
freq_table(bsl, c("status", "sample"))
freq_table(bsl, "genes")

table(bsl[, c("status", "sample")])
glist <- levels(as.factor(bsl$genes))
#' The target gene is RPL23.
#'
#' ### Normalize against reference gene
#' Pros: no need for accurate quantification of starting material
#' Cons: require reference gene(s) with stable expression levels
#'
#' | |Target gene | Reference gene|
#' | --- | --- | --- |
#' |Callibrator|$C_T(target, calibrator)$|$C_T(reference, calibrator)$|
#' |Test|$C_T(target, test)$|$C_T(reference, test)$|
#'
#' The most common approaches are:
#'
#' - Livak, aka -$\Delta \Delta C_T$ method
#' - $\Delta C_T$ method using a reference gene
#' - Pfaffl method
#'
#' ## Delta-delta method $2^{-\Delta \Delta Ct}$
#' Steps:
#'
#' 1. Normalize $C_T$(target gene) to $C_T$(reference gene)
#' 2. Noramlize $\Delta C_T$ of test sample to $\Delta C_T$ of calibrator
#' 3. Calculate expression ratio, or fold change
#'
#' Using the example dataset:
#'
#' 1. Taking the mean of the multiple measures
m_ct <- bsl %>%
  group_by(ID, status, sample, genes) %>%
  dplyr::summarise_at(dplyr::vars(ct), ~ mean(.x, na.rm = T))

head(m_ct)

#' 2. calculate  the first $\Delta Ct$ (normalize against reference gene) in
#' each sample.
#' $\Delta C_T = -(C_T(reference) - C_T(target) )= C_T(target) - C_T(reference)$
# pivot long to wide table for calculation
d_ct_sample <- m_ct %>%
  tidyr::pivot_wider(
    id_cols = c("ID", "status", "sample"),
    names_from = "genes", values_from = "ct"
  ) %>%
  dplyr::arrange(ID, status, sample) %>%
  mutate_at(glist, ~ purrr::map_dbl(., function(x) {
    x - RPL23
  })) %>%
  select(-RPL23) %>%
    data.frame()

# plot value distribution
glist2 <- make.names(glist[-12]) # remove reference gene
plot(density(d_ct_sample[, glist2[1]]),
  ylim = c(0, 0.8),
  xlim = c(0, 20),
  xlab = "", main = expression(paste(Delta, C[T], "  distributions"))
)
invisible(lapply(seq(2, length(glist2)), function(i) {
  g <- glist2[i]
  lines(density(d_ct_sample[, g], na.rm = T),
    xlab = "", ylab = "", main = "",
    col = i, lty = i,
    add = T
  )
}))
legend("topleft", glist2,
  lty = 1:14,
  col = 1:14, cex = 0.5, bty = "n"
)

#' 3. Normalize $\Delta C_T$ of disease sample to $\Delta C_T$  of Healthy samples
#' $\Delta C_{T, D} - \Delta C_{T, H}$ on the populational level. So I take
#' the mean of the expression levels in diseases and healthy subjects.
diff <- function(x) x - first(x) # first row is healthy subjects
dd_ct <- d_ct_sample %>%
  dplyr::arrange(status) %>%
  group_by(status, sample) %>%
  summarise_at(glist2, ~ mean(.x, na.rm = T)) %>%
  # healthy subject on the top row
  dplyr::arrange(desc(status)) %>%
  ungroup() %>%
  mutate_at(glist2, list("chg" = diff))

dd_ct_glist <- grep("chg", names(dd_ct), value = T)

dd_ct %>%
  mutate_at(dd_ct_glist, ~ purrr::map_dbl(., function(x) 2^(-x))) %>%
  select(status, sample, ends_with("chg")) %>%
  t() %>%
    DT::datatable()

#' ## $\Delta C_T$ method
#' $\Delta C_T$ method calcluates the difference of gene expression
#'  between the reference and target
#' $C_T$ values for each sample. It uses an expression value of the calibrator
#' sample that is not 1.
#' Steps:
#'
#' 1. Normalize target gene for each sample with
#' $2^{C_T(reference)- C_T(target)} = 2^{- \Delta C_T}$
#' 2. Determine ratio of expression: Control expression = control/ control;
#' Disease expression = disease/control
#' Use the example data, the first step was already half done
exp <- d_ct_sample %>%
  mutate_at(glist2, ~ purrr::map_dbl(., function(x) 2^(-x)))

plot(density(exp[, glist2[1]]),
  ylim = c(0, 1000),
  xlim = c(-0.06, 0.25),
  xlab = "",
  main = expression(paste("Gene expression level distributions " * 2^-Delta^C[T]))
)
invisible(lapply(seq(2, length(glist2)), function(i) {
  g <- glist2[i]
  lines(density(exp[, g], na.rm = T),
    xlab = "", ylab = "", main = "",
    col = i, lty = i,
    add = T
  )
}))
legend("topright", glist2,
  lty = 1:14,
  col = 1:14, cex = 0.5, bty = "n"
)

#' ### Summary statistics on $2^{-\Delta C_T}$ - relative gene expression level in each group
#' The gene expression levels are relative to the gene expression level to the reference gene.
exp %>%
  group_by(status, sample) %>%
  get_summary_stats(type = "mean_se") %>%
  mutate(n_mean_se = paste0("N = ", n, ", ", mean, " (", se, ")")) %>%
    select(-c(n, mean, se)) %>%
    tidyr::pivot_wider(
      names_from = c("status", "sample"),
      values_from = c(n_mean_se)
    ) %>%
    DT::datatable(rownames = F)


#' ### Comparison between two groups
#' The comparison between two groups/conditions or against referece group
#' is done by taking the ratio of the expression levels, i.e.
#' $$\frac{exp_{disease}}{ exp_{ref}} =\frac{2^{\Delta C_T(D)} }{2^{\Delta C_T (reference)}}$$.
#'
#' This is enssentially the same as the delta delta method.
#'
#' ## Pfaffl method
#' This is used when reaction efficiencies of the target and reference genes are
#' not similar.
#' $$Ratio = \frac{E_{target}^{[\Delta C_T, target(calibrator - test)]}}{E_{ref}^{[\Delta C_T, ref(calibrator - test)]}}$$
#'
#' # Analysis using R package

#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/061_qPCR.R', output_dir = 'output')
# knitr::spin('src/061_qPCR.R', format = 'Rmd', knit = FALSE)