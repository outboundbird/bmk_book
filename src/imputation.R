#' ---
#' title: Simple illustration for detection limit problem
#' subtitle: 'SAR: NA , Study: simulation study'
#' author:  Siying Huang (E0482362), Biomarker statistics team
#' date: 'created: 2023-09-07 , updated (`r Sys.Date()`)'
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
# knitr::opts_knit$set(root_dir='/mnt/c/Users/e0482362/Work/bookdown/bookdown/src')
knitr::opts_chunk$set(echo = T, comment = "", message = F, warning = F, error = F)
options(width = 100)
#+ libs
library(here)
library(NADA)
library(ggplot2)
library(dplyr)
# create groups
# y ~ b *g + e
# detection limit at 20% at 50%
# method of handling detection limit
# bias in effect estimate, type 2 error
#' # Simulation theme
#' In this simulation of the detection limit scenario, I generate data with following parameters:
#'
#' - Two data sets of comparison groups, one balanced with ratio of 1:1, the other with unbalanced group with ratio of 1:3. There are total of 40 subjects in each dataset.
#' - missing percentage. For both sets, the detection limit are set to be the 25% and 50% of the distribution on Y.
g1 <- rep(c(1,2), 20) %>% factor(labels = c('A','B'))
g2 <- c(rep(1, 10), rep(2, 30)) %>% factor(labels = c('C','D'))
set.seed(1234)
b0 <- 30
b1 <- 3
y1 <- b0 + b1 * as.integer(g1) + rnorm(40, 5, 2)
y2 <- b0 + b1 * as.integer(g2) + rnorm(40, 5, 2)

summary(y1)
summary(y2)

# set up missing values at 20% and 50%
dl1 <- quantile(y1, c(0.2, 0.5))
dl2 <- quantile(y2, c(0.2, 0.5))

y1_20 <- ifelse(y1 <= dl1[1], NA, y1)
y1_50 <- ifelse(y1 <= dl1[2], NA, y1)

y2_20 <- ifelse(y2 <= dl2[1], NA, y2)
y2_50 <- ifelse(y2 <= dl2[2], NA, y2)

#' # Missing pattern
#' ## Dataset 1
table(is.na(y1_20), g1)
table(is.na(y1_50), g1)

#' ## Dataset 2
table(is.na(y2_20), g2)
table(is.na(y2_50), g2)

#' # Effect estimate with the presence of detection limit {.tabset}
#' ## Dataset1
par(mfrow = c(1, 3))
boxplot(y1 ~ g1, ylim = c(30, 50), main = "Full data")
boxplot(y1_20 ~ g1, ylim = c(30, 50), main = "25% under DL")
abline(h = dl1[1], lty = 2, col = 2)
boxplot(y1_50 ~ g1, ylim = c(30, 50), main = "50% under DL")
abline(h = dl1[2], lty = 3, col = 2)
mtext("Dataset 1", side = 3, line = -1, outer = T)

rbind(
  summary(lm(y1 ~ as.factor(g1)))$coeff[2, ],
  summary(lm(y1_20 ~ as.factor(g1)))$coeff[2, ],
  summary(lm(y1_50 ~ as.factor(g1)))$coeff[2, ]
) %>%
  DT::datatable(caption = "Dataset1") %>%
  DT::formatRound(c(1, 2, 3)) %>%
  DT::formatSignif(4)

#' ## Dataset2
boxplot(y2 ~ g2, ylim = c(30, 50), main = "Full data")
boxplot(y2_20 ~ g2, ylim = c(30, 50), main = "25% under DL")
abline(h = dl2[1], lty = 2, col = 3)
boxplot(y2_50 ~ g2, ylim = c(30, 50), main = "50% under DL")
abline(h = dl2[1], lty = 2, col = 3)
mtext("Dataset 2", side = 3, line = -1, outer = T)

rbind(
  summary(lm(y2 ~ as.factor(g2)))$coeff[2, ],
  summary(lm(y2_20 ~ as.factor(g2)))$coeff[2, ],
  summary(lm(y2_50 ~ as.factor(g2)))$coeff[2, ]
) %>%
  DT::datatable(caption = "Dataset2") %>%
  DT::formatRound(c(1, 2, 3)) %>%
  DT::formatSignif(4)



#' # Methods for detection limit problem
#' Here we compare two simple methods: replacement at detection limit and ROS
hist(y1, breaks = length(y1))
hist(y2, breaks = length(y2))

ds1 <- list(
  # full = data.frame(y = y1, g = g1),
  y20 = data.frame(y = y1_20, g = g1, ind= is.na(y1_20)),
  y50 = data.frame(y = y1_50, g = g1, ind = is.na(y1_50))
)

ds2 <- list(
  # full = data.frame(y = y2, g = g2),
  y20 = data.frame(y = y2_20, g = g2, ind = is.na(y2_20)),
  y50 = data.frame(y = y2_50, g = g2, ind = is.na(y2_50))
)

# replacement

#
cenboxplot(ds1$y20$y, ds1$y20$ind, ds1$y20$g)
cendiff(ds1$y20$y, ds1$y20$ind, ds1$y20$g)
#' <details><summary>Session Info</summary>
sessionInfo()
#' </details>
#+ echo = F, eval = F
# Markdown --------------------------------------------------------
# rmarkdown::render('src/imputation.R', output_dir = 'output')
# knitr::spin('src/imputation.R', format = 'Rmd', knit = FALSE)