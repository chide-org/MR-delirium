outcome <- fread("D:/A-MR/outcome/rs转换_delirium.csv")

outcome$id <- "delirium"
outcome$phe <- "delirium"
write.csv(outcome, file = "D:/A-MR/outcome/complete_delirium.csv")

read_outcome <- read_outcome_data(
    "D:/A-MR/outcome/complete_delirium.csv",
    snps = NULL,
    sep = ",",
    phenotype_col = "phe",
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    eaf_col = "af_alt",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    pval_col = "pval",
    gene_col = "nearest_genes",
    id_col = "id",
    log_pval = T,
    chr_col = "#chrom",
    pos_col = "pos"
)
viewn(read_outcome)
fwrite(read_outcome, file = "outcome.csv", nThread = 8)
# quote为T时，字符型数据会被引号包裹，na为NA时，缺失值会被写为NA，qmethod为double时，双引号将会被双引号包裹，nThread为线程数
# nThread不同的值时，写入相同文件的时间（我的cpu为物理4核，逻辑8核）
# 3:10.01s
# 4:11.65s
# 5:6.17s
# 6:6.48s
# 7:7.95s
# 8:5.09s
# 9:5.22s
# 10:6.49s
# 11:6.17s
