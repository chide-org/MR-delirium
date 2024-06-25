library(TwoSampleMR)
# outcome
outcome <- fread("outcome/zhanwang")

library(dplyr)
oucome_empty <- outcome %>%
  filter(rsids == "") %>%
  select("#chrom", "pos", "ref", "alt", "rsids")
fwrite(oucome_empty, "outcome/oucome_empty.csv")
outcome_rename <- outcome %>%
  select(SNP, Pvalue, Effect_allele, Non_Effect_allele, Beta, SE, EAF) %>%
  mutate(id = "delirium")

write.csv(outcome_rename, file = "outcome/delirium.csv")
df <- outcome[1:23, c(SNP, Pvalue, Effect_allele)]
viewn(df)
# ieu
ieus <- c("prot-a-530")
ieu_id <- c("ieu-b-110")
View(exposure)
# for
for (ieu_id in ieus) {
  # exposure
  exposure <- extract_instruments(outcomes = ieu_id)
  dat <- merge(exposure, outcome, by.x = "SNP", by.y = "MarkerName")
  write.csv(dat, file = "temp.csv")
  dat_read <- read_outcome_data(
    snps = exposure$SNP,
    filename = "temp.csv",
    sep = ",",
    snp_col = "SNP",
    pval_col = "Pvalue",
    effect_allele_col = "Effect_allele",
    other_allele_col = "Non_Effect_allele",
    beta_col = "Beta",
    se_col = "SE",
    eaf_col = "EAF",
  )
  file.remove("temp.csv")
  # harmon
  dat_harmon <- harmonise_data(
    exposure_dat = exposure,
    outcome_dat = dat_read
  )
  # mr
  mr <- mr(dat_harmon)
  all_mr <- rbind(all_mr, empty_row, mr)
}
write.csv(all_mr, file = "ADoutcome.csv")


# init
empty_row <- data.frame(matrix(NA, nrow = 1, ncol = ncol(mr)))
colnames(empty_row) <- colnames(mr)
all_mr <- mr
all_mr <- all_mr[1, ]
