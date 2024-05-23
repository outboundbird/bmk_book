library(dplyr)
# qPCR -----------------------------------------
data_long <- haven::read_sas(file.path("data", "pf_h.sas7bdat")) %>%
  data.frame() %>%
  select(USUBJID, PFSTRESN, VISIT, PFRESCAT, PFGENRI, PFREFID, PFSPID) %>%
  arrange(USUBJID, VISIT, PFGENRI)
head(data_long, 10)
dim(data_long)
data <- readRDS("data/bmk-qtpcr_processed.Rds")
id_list <- setNames(unique(data_long$USUBJID), seq(41))
switch_id <- function(id, id_list) {
  # id_list must have names for each element in the list
  if (id %in% id_list) {
    names(id_list[id_list == id])
  }
}
data %>%
  select(USUBJID, DSGROUP, AGE, SEX) %>%
  distinct() %>%
  right_join(data_long, by = c("USUBJID")) %>%
  mutate(
    PFSTRESN = PFSTRESN + rnorm(1),
    status = if_else(DSGROUP == "HEALTHY SUBJECTS", "Healthy", "Disease"),
    sample = case_when(
      PFRESCAT == "LESIONAL NON-PRURITIC SKIN" ~ "Lesion N",
      PFRESCAT == "LESIONAL PRURITIC SKIN" ~ "Lesion P",
      PFRESCAT == "NON-LESIONAL SKIN" ~ "Normal",
      PFRESCAT == "POST-LESIONAL HEALED SKIN NON-PRURITIC NON LESIONAL" ~ "Post_lesion"
    ),
    ct = PFSTRESN,
    pos = PFREFID,
    sponsor = PFSPID,
    genes = gsub(" GENE", "", PFGENRI),
    Day = VISIT,
    ID = sapply(USUBJID, function(x) switch_id(x, id_list))
  ) %>%
  select(ID, Day, status, sample, pos, sponsor, genes, ct, AGE, SEX) %>%
  saveRDS("data/qPCR.rds")
