library(dplyr)
library(ggpubr)

visit <- rep(c("baseline", "day15", "day29"), each = 6)
base_norm <- rnorm(6, 1)
base_lesion <- rnorm(6, 10)
day15_imp <- rnorm(6, 7)
day15_plc  <- rnorm(6, 10)
day29_imp <- rnorm(6, 4)
day29_plc <- rnorm(6, 9.5)

df <- data.frame(
  base_norm,
  base_lesion, day15_imp, day15_plc, day29_imp, day29_plc
) %>%
tidyr::pivot_longer(1:6, names_to = "time", values_to = "geneExp") %>%
  mutate(visit = stringr::str_extract(time, "\\w+(?=\\_)"),
  group = stringr::str_extract(time, "(?<=\\_)\\w+"))

ggline(df,
  x = "visit",
  y = "geneExp",
  group = "group",
  color = "group",
 add = "mean_se",
 error.plot = "pointrange"
)
