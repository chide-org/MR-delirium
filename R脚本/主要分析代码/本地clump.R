# 本地去除连锁不平衡的SNP
gc()
library(plinkbinr)
library(ieugwasr)

clump <- "D:/A-MR/clump_tool"
setwd(clump)
getwd()
# 教程：https://mp.weixin.qq.com/s/TDEZYlp2AHHA8Yqld8GKHA


# 输入待clump的数据
data_clumped <- ld_clump(
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 0.99,
  pop = "EUR",
  dplyr::tibble(rsid = dat$SNP, pval = dat$p, id = dat$id),
  plink_bin = "plink_win64_20231211/plink.exe",
  bfile = "1kg.v3/EUR"
)
viewn(data_clumped)
# 转换并合并数据
dat <- data_clumped %>%
  rename(SNP = rsid, pval.exposure = pval, id.exposure = id) %>%
  left_join(dat, by = c("SNP", "pval.exposure", "id.exposure"))

