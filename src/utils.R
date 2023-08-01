denote_p <- function(p) {
  thresh <- setNames(
    c(1, 0.1, 0.05, 0.01, 0.001, 0.0001),
    c("ns", ".", "*", "**", "***", "****")
  )
  names(thresh)[max(which(p < thresh))]
}


power <- function(base) {
  function(log) {
    round(base^log, 2)
  }
}
exp2 <- power(2)


t_test <- function(data, y, contrast, subset_var, ...) {
  formu <- as.formula(paste(y, "~", contrast))
  sub_idx <- data[, subset_var]
  contr_grp <- data.frame(table(data[, c(subset_var, contrast)])) %>%
    filter(.data[[subset_var]] == T & Freq > 0) %>%
    select(.data[[contrast]])

  contr_grp <- paste(as.character(contr_grp[, 1]), collapse = ": ")

  rst <- t.test(formu, data,
    subset = sub_idx,
    na.action = na.omit,
    ...
  )

  diff <- round(rst$estimate[1] - rst$estimate[2], 2)
  lcl <- round(rst$conf.int[1], 2)
  ucl <- round(rst$conf.int[2], 2)
  p <- round(rst$p.value, 4)
  se <- round(rst$stderr, 2)
  sig <- denote_p(p)
  setNames(
    c(y, contr_grp, diff, se, lcl, ucl, p, sig),
    c("var", "contrast", "diff", "se","lcl", "ucl", "p.val", "p_sig")
  ) %>%
    t()
}
