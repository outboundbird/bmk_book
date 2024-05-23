#' ---
#' title: Sample size calculation for RNA sequencing analysis
#' subtitle: 'SAR: SAR440340  , Study: NA'
#' author:  Siying Huang (E0482362), Biomarker statistics team
#' date: 'created: 2023-12-11 , updated (`r Sys.Date()`)'
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
#+ libs
library(dplyr)
library(kableExtra)
library(ssizeRNA)
library(RnaSeqSampleSize)
library(RNASeqPower)
# https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS4265
#' # Experiment setup
#'
#' - **Variations.** Assuming biological coefficient of variation (CV) in comparison groups are the same, for human beings, the biological CV ranges from 0.32 to 0.74 witha median of 0.43. The CV is also tissue type dependent. In this exercise the *CV* is set to **0.4** assuming the variations among the comparison groups are the same. This is relatively conservative because the genes from lung tissue are much less mapped (due to tissue-specific expression) and biological CV is relatively high. [See this article for details](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3842884/) Another type of variation comes from the machine, this can reduced to some extent by increasing sequencing depth. In this exercise, the *sequencing depth* is set to be **10**.
#'
#' - **Dispersion** Dispersion in the variance for negative binomial distribution. In this exercise, the dispersion factors are set to be **0.1** and **0.2**.
#'
#' - **Significance threshold.** Taken the multiple comparison issue into consideration for significance threshold definition and assuming after filtering, there will be 15,000 to 18,000 genes left for analysis, the most conservative p-value is calculated as $\alpha=$ 0.05/$N_{tests} =$ 3.3e-06.Note that the number of gene to be tested is derived from PBMC samples. The number varies depending on the tissue type. Some genes may are expressed in certain tissues. In this exercise the number of genes to be tested is set to 15,000.
#'
#' - **Fold change.** The ratio of gene expression levels in two different status/ groups. The fold changes in exercise are set to be **1.5, 1.75** and **2**, which correspond to 50%, 75% and 100% of up-regulation.
#'
#' - **Statistical power.** The probability to reject a gene is not differentially expressed when the gene is truly differentially expressed. By default we want to have at least **80%** probabilty to make the right call.
#'
#' - **Fraction of non-differentially expressed genes.** The fraction of genes out of total gene pool that are differentially expressed between the comparison groups/status. e.g. 10% of 18,000 genes would be differentially expressed results in 1800 genes. In this exercise, this fraction is set to be 20%. This is completely a guess.
#'


fc <- c(1.5, 1.75, 2)
n_gene <- 15000
psig <- 0.05 / n_gene
fdr_thresh <- 0.05
seq_depth <- c(10, 20)

#' # Sample size
#' ## Simple sample size calculation
#' The sample size was calculated with most conservative false positive control.
#' This method considered CV, sequencing depth as factors account for the variations in the observations. The CVs are set to 0.4, 0.5 in columns, sequencing depth to be 10 and 20 in rows. The fold changes are expected to be 1.5, 1.75 and 2 respectly.

tab <- RNASeqPower::rnapower(
  depth = c(10, 20), cv = c(0.4, 0.5),
  effect = fc, alpha = psig, power = 0.8
) %>%
  data.frame(check.names = F, fix.empty.names = F)

header <- names(tab)
cv <- paste0("CV=", stringr::str_sub(header, 1, 3))
# fc <- stringr::str_sub(header, 5, 8)
names(tab) <- cv

kable(tab, digit = 0) %>%
  add_header_above(c(" " = 1, "FC=1.5" = 2, "FC=1.75" = 2, "FC=2" = 2))

#' ## Power estimation from proposed sample size
#'
#' The sample sizes derived from the above are used to estimate statistical
#' powers are using simulation assuming
#' 80% of genes are not differentially expressed, and the dispersion factors are
#' 0.1 , 0.2 and 0.5. This simulaiton does not directly consider the variations from
#' biological or technical duplicates.
#'
#'
#+ fc1.5, cache = T
dp <- c(0.1, 0.2, 0.5, 0.1, 0.2, 0.5, 0.1, 0.2, 0.5)
n <- c(26, 26, 26, 50, 50, 50, 128, 128, 128)
sim <- mapply(function(m, dp) {
  rst <- ssizeRNA::check.power(15000,
    pi0 = 0.8,
    m = m, mu = 10,
    disp = dp, fc = 1.5, sims = 10
  )
  c(rst$pow_bh_ave, rst$pow_qvalue_ave, rst$fdr_bh_ave, rst$fdr_qvalue_ave)
}, n, dp,
SIMPLIFY = T, USE.NAMES = T
)

sim <- data.frame(sim)

names(sim) <- dp
rownames(sim) <- c("power_BH", "power_ST", "fdr", "qvalue")

kable(sim, digits = 2, caption = 'FC = 1.5') %>%
  add_header_above(c(" ", "N =26" = 3, "N=50" = 3, "N=128" = 3))

#+ fc2, cache = T
sim <- mapply(function(m, dp) {
  rst <- ssizeRNA::check.power(15000,
    pi0 = 0.8,
    m = m, mu = 10,
    disp = dp, fc = 2, sims = 10
  )
  c(rst$pow_bh_ave, rst$pow_qvalue_ave, rst$fdr_bh_ave, rst$fdr_qvalue_ave)
}, n, dp,
SIMPLIFY = T, USE.NAMES = T
) %>%
data.frame(check.names = F)

names(sim) <- dp
rownames(sim) <- c("power_BH", "power_ST", "fdr", "qvalue")

kable(sim, digits = 2, caption = 'FC = 2') %>%
  add_header_above(c(" ", "N =26" = 3, "N=50" = 3, "N=128" = 3))

#+ echo = F
# RnaSeqSampleSize::

#+ data_source, echo = F
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212331
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175829

#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/samplesize.R', output_dir = 'output')
# knitr::spin('src/samplesize.R', format = 'Rmd', knit = FALSE)
