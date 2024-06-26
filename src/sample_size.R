#' ---
#' title: "Sample size calculation for nanostring validation study"
#' author: "Siying Huang (E0482362), Vinh Truong (E0437543)"
#' date: "Jan 19 2022 update (`r Sys.Date()`)"
#' output:
#'   html_document:
#'     css: !expr here::here("notebooks/sanofi_doc.css")
#'     code_folding: "hide"
#'     toc: yes
#'     toc_float:
#'       collapse: no
#' ---
#'
#+ setup, include = FALSE
knitr::opts_knit$set(
  root.dir = here::here(), #"/mnt/c/Users/e0482362/Work/cd40l/CD40L_transcriptomic/",
  comment = "",
  message = F,
  warning = F
)

library(ICC.Sample.Size)
library(tidyr)
library(kableExtra)

#' # Data structure example
#+
sample <- rep(c("RNAseq", "nanostring"), 3)
genes <- matrix(rnorm(30), ncol = 6)
rownames(genes) <- paste0("gene", c(1:5))
colnames(genes) <- sample
kable(genes) %>%
  add_header_above(c("", "ID1" = 2, "ID2" = 2, "ID3" = 2))

#' ## Objectives of validation study
#' - assess the agreement of two measurements
#' - construct a calibration equation in order to unify the cut off point between two measures.
#'
#' # Agreement between nanostring and RNASeq measurements
#' ## Inter-class correlation coefficient
#'
#' - Pearson correlation
#' - Spearman correlation
#'
#' $$ \frac{\sum_i\rho(G_{j}(Rseq), G_{j}(nano)))}{n_i}
#' $$
#' - i - # of genes
#' - j - # of subjects
#'
#' Need to visually examine if there exist ceiling and flooring effect
#' between two technologies
#'
#' ## Intraclass Correlation Coefficient (ICC)
#' Generally speaking, the ICC determines the reliability of ratings by
#' comparing the variability of different ratings of the same individuals
#' to the total variation across all ratings and all individuals.
#'
#' - A high ICC (close to 1) indicates high similarity between values from
#' the same group.
#' - A low ICC (ICC close to zero) means that values from the
#' same group are not similar.
#'
#' The sample size is calculated based on following parameter:
#'
#' - ICC0 - ICC of null hypothesis (e.g. 0 no agreement or a value based on previous observations)
#' - ICC - hypothetical ICC
#' - alpha at 0.05
#' - k = 2 (number of raters: nanostring vs. RNAseq)
#' - two tails
#'
#' ### Reference
#' [1]T. K. Koo and M. Y. Li, “A Guideline of Selecting and Reporting Intraclass Correlation Coefficients for Reliability Research”, Journal of Chiropractic Medicine, vol. 15, no. 2, pp. 155–163, Jun. 2016, doi: 10.1016/j.jcm.2016.02.012.
#'
#' [2]S. D. Walter, M. Eliasziw, and A. Donner, “Sample size and optimal designs for reliability studies”, Statistics in Medicine, vol. 17, no. 1, pp. 101–110, 1998, doi: 10.1002/(SICI)1097-0258(19980115)17:1<101::AID-SIM727>3.0.CO;2-E.
#'
#' ### Power at 80%, single tail
k2a005p08 <- calculateIccSampleSize(by = "both", tails = 1)
table_n <- k2a005p08[[2]]

knitr::kable(table_n,
caption = "Table of sample size (N) given null and alterative hythosis of ICC.") %>%
  column_spec(2:12, color = ifelse(rownames(table_n) > 0.5, "red", 1)) %>%
  add_header_above(c("", "ICC0" = 20)) %>%
  pack_rows("ICC", 1, 20)

#' ###  Power at 90%, single tail
#'
k2a005p09 <- calculateIccSampleSize(power = 0.9, by = "both", tails = 1)
table_n <- k2a005p09[[2]]

knitr::kable(table_n,
  caption = "Table of sample size (N) given null and alterative hythosis of ICC."
) %>%
  column_spec(2:12, color = ifelse(rownames(table_n) > 0.5, "red", 1)) %>%
    add_header_above(c("", "ICC0" = 20)) %>%
    pack_rows("ICC", 1, 20)

#' # Considerations for calibration of two measures
#'
#' ## Reference:
#' [3] L. Tian, R. A. Durazo-Arvizu, G. Myers, S. Brooks, K. Sarafin, and C. T. Sempos, “The estimation of calibration equations for variables with heteroscedastic measurement errors,” Statistics in Medicine, vol. 33, no. 25, pp. 4420–4436, 2014, doi: 10.1002/sim.6235.
#'
#' ## Calibration function
#' The goal is to associate a test value (e.g. RNAseq) to a reference value (e.g. nanostring) via a calibration function.
#' We collect N samples of paired measured:
#' $\{(X_i, Y_i), i = 1... N\}$, where
#' $$ X_i = \tilde{X_i} + \epsilon_i^x \\
#' Y_i = \tilde{Y_i} + \epsilon_i^y\\
#' \tilde{Y_i} = f_0(\tilde{X_i}) +e_i
#' $$
#' where $\epsilon_i^x ,\epsilon_i^y, e_i$ are zero mean errors with variance
#' $\sigma_i^{x2}, \sigma_i^{y2}, \sigma_0^2$
#'
#' In simple calibration function, we assume a linear function: $f_0(x) = \beta_0+\beta_1X$.
#' We want to estimage $\beta_0, \beta_1$ to solve the calibration function.
#'
#' Some considerations for calibration function:
#'
#'  - the variance of the measurement error often depends on the underlying level.
#' - standard deviation for the measurement error is proportional to the underlying value
#' - the CV (SD/mean) of themeasurement error is **approximately a constant**
#'
#' therefore the estimate of variance, assuming:
#' $$ \sigma_i^x = CV_x \times \tilde{X} \\
#' \sigma_i^y = CV_x \times \tilde{Y} \\
#' $$
#'
#' For estimating parameters in calibration function see the reference paper for different senarios.
#'
#' ## Sample size
#'
#' $$
#' N = max_{k=1..K}\bigg\{\frac{(2Z_{1-\alpha})^2\times \sigma^2(x_k)}{\delta_k^2}\bigg\}\\
#' \sigma^2(x_k) =(1, x_k)\sum\bigg(1, x_k\bigg)
#' $$
#'
#' $\sum$ is variance, co-variance matrix for estimated intercept and slope of the calibration function.
#' $\delta_k, k = 1...K$ precision levels.
#'
#' Parameters to consider before calculation:
#'
#' - x0: x value plan to calibrate with estimated calibration equation (e.g. $\beta_0=0, \beta_1=1$, and this never happens in reality)
#' - d0: 95% CI of calibrated x value
#' - x:  emperical observation of targeted distribution (e.g. mRNA levels from RNA seq or from nanostring)
#' - CVx: coefficient variation of measurement X (e.g. nanostring), assume to be constant
#' - CVy: coefficient variation of measurement Y (e.g. RNAseq), assume to be constant
#'
#' ## Sample size calculation {.tabset}
#' ### CVx=CVy = 0.1
#'
#+ fig.dim = c(9, 7), results ='hide'
x0 <- c(40, 50, 60, 70, 80, 90)
d0 <- c(5, 10, 20, 30, 40, 100)

calc_samplesize <- function(x0, d0, x =seq(5, 1000, length = 1000), CVx =0.1, CVy=0.1) {
  CVcalibration::samplesize(x0, d0,
    x = x, # assuming RNAseq raw count data
    CVx = CVx, CVy = CVy
  )
}
par(mfrow = c(2, 3))
lapply(d0, function(d) {
  rst <- sapply(x0, function(x) calc_samplesize(x0 = x, d0 = d))
  tab <- do.call("rbind", rst)
  # rownames(tab) <- x0
  ymin <- min(tab)
  ymax <- max(tab)
  plot(x0, tab[, 1],
    type = "b",
    ylim = c(ymin, ymax),
    ylab = "N",
    main = paste("d0=", d)
  )
  lines(x0, tab[, 3], lty = 3, type = "b", pch = 2)
  legend("topright", c("1 CV known", "both CVs known"),
    bty = "n",
    lty = c(1, 3), pch = c(1, 2)
  )
  return(NULL)
})

#' ### CVx=0.2, CVy=0.2
#'
#+ fig.dim = c(9, 7), results ='hide'

par(mfrow = c(2, 3))
lapply(d0, function(d) {
  rst <- sapply(x0, function(x) calc_samplesize(x0 = x, d0 = d, CVx=0.2, CVy = 0.2))
  tab <- do.call("rbind", rst)
  # rownames(tab) <- x0
  ymin <- min(tab)
  ymax <- max(tab)
  plot(x0, tab[, 1],
    type = "b",
    ylim = c(ymin, ymax),
    ylab = "N",
    main = paste("d0=", d)
  )
  lines(x0, tab[, 3], lty = 3, type = "b", pch = 2)
  legend("topright", c("1 CV known", "both CVs known"),
    bty = "n",
    lty = c(1, 3), pch = c(1, 2)
  )
  return(NULL)
})


#' ### CVx=0.1, CVy=0.05
#'
#+ fig.dim = c(9, 7), results ='hide'

par(mfrow = c(2, 3))
lapply(d0, function(d) {
  rst <- sapply(x0, function(x) calc_samplesize(x0 = x, d0 = d, CVy = 0.05))
  tab <- do.call("rbind", rst)
  # rownames(tab) <- x0
  ymin <- min(tab)
  ymax <- max(tab)
  plot(x0, tab[, 1],
    type = "b",
    ylim = c(ymin, ymax),
    ylab = "N",
    main = paste("d0=", d)
  )
  lines(x0, tab[, 3], lty = 3, type = "b", pch = 2)
  legend("topright", c("1 CV known", "both CVs known"),
    bty = "n",
    lty = c(1, 3), pch = c(1, 2)
  )
  return(NULL)
})


#' ### CVx=0.05, CVy=0.05
#'
#+ fig.dim = c(9, 7), results ='hide'

par(mfrow = c(2, 3))
lapply(d0, function(d) {
  rst <- sapply(x0, function(x) calc_samplesize(x0 = x, d0 = d, CVx = 0.05, CVy = 0.05))
  tab <- do.call("rbind", rst)
  # rownames(tab) <- x0
  ymin <- min(tab)
  ymax <- max(tab)
  plot(x0, tab[, 1],
    type = "b",
    ylim = c(ymin, ymax),
    ylab = "N",
    main = paste("d0=", d)
  )
  lines(x0, tab[, 3], lty = 3, type = "b", pch = 2)
  legend("topright", c("1 CV known", "both CVs known"),
    bty = "n",
    lty = c(1, 3), pch = c(1, 2)
  )
  return(NULL)
})

#' ## Choosing samples range
#' > There are two important considerations for choosing sampling methods: (i) the
#' > selected Xi values should cover the entire region of interest and (ii) the particular sampling scheme should
#' > help us to obtain an accurate estimate for the calibration function. There are several obvious choices:
#'
#' > (1) randomly sampling from the stored samples (sample more from important regaions, but less accurate)
#' >
#' > (2) uniformly sampling from the given interval of interest (more accurate, less samples from important regions to investigate linear assumption) ;
#' >
#' > (3) a hybrid of the aforementioned two sampling methods (balance between above two):
#' >
#' > (i) dividing the range of interest into intervals using the quantiles of the observed samples
#' > (ii) uniformly sampling equal number of Xi’s  from each of the subintervals.
#'
#+ fig.dim = c(8,3)
par(mfrow = c(1, 3))
set.seed(789)
x <- rnorm(500, 2)
x1 <- runif(500, -1, 5)
qt <- quantile(x)
x2 <- c(runif(100,qt[1], qt[2]), runif(100, qt[2], qt[3]), runif(100, qt[3], qt[4]), runif(100,qt[4], qt[5]))
hist(x, breaks =30, main ="random sampling")
hist(x1, breaks = 30, main ="uniform sampling")
hist(x2, breaks = 30, main = "hybrid sampling")

#' ## Estimating coefficient variation (CV) within the same measures
#' - repeated measures in the same sample to obtain the estimation of CV in both techonologies
#'
#'
#'
#' ## Consideration for prediction precision
#' Calibration can be viewed as a prediction model. Without any error measurement and assuming
#' constant error variance, a linear model can be used:
#' $$Y_i = a\beta X_i +b + \epsilon$$
#'
#' Uncertainy will be considered when translating the cut-off from RNAseq data
#' to Nanostring. It boils down to estimate with prediction the residual variance.
#'  When a sample of size $n$ is drawn from a normal
#' distribution, a $1-\alpha$ CI of the unknown population variance is given by
#' $$\frac{n-1}{\chi^2_{1-\alpha/2, n-1}}s^2  < \sigma^2 < \frac{n-1}{\chi^2_{\alpha/2, n-1}}s^2$$
#' The fold change between the bounds (multiplicative margin of error) can be
#'  used as a precision proxy [4].
#'
#' Formula can be adapted with $p$ parameters [5]. Non linear relationship
#' can be considered with a restricted cubic splines and k knots (k-1 parameters)
#' We assume that k = 3



dat <- data.frame(
  MMOE = c(1.1, 1.2),
  linear = c(235, 71),
  non_linear = c(238, 74)

)

knitr::kable(dat)
#' ### Reference
#' [4]F. E. Harrell Jr., “Regression Modeling Strategies”. Springer International Publishing, 2016
#'
#' [5]R. D. Riley, J. Ensor, K. I. E. Snell, et al. “Calculating the sample size required for developing a clinical prediction model”. BMJ. 2020

#'
#' # Final consideration for sample size
#' The sample size consideration should serve the purpose of analysis, which in our case is calibraion
#' and agreement test. These two objectives might render different numbers. In order to gurantee
#' the quality of prelminary analysis, we should take the max{N1, N2, N3}.
#'



# rmarkdown::render("src/sample_size/sample_size.R")
