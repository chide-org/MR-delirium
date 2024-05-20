# 本地去除连锁不平衡的SNP
library(plinkbinr)
library(ieugwasr)
library(dplyr)
library(data.table)


# 本地clump参照https://mp.weixin.qq.com/s/TDEZYlp2AHHA8Yqld8GKHA



# 读取文件
data <- fread("exposure/Evangelou2018 DBP.txt")
nrow(data)
viewn(data)
# 过滤数据
data_subset <- subset(data, P < 5e-8)
nrow(data_subset)
colnames(data_subset)
viewn(data_subset)
data_subset$id <- "id"
# 输入待clump的数据
data_clumped <- ld_clump(
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 0.99,
  pop = "EUR",
  dplyr::tibble(rsid = data_subset$MarkerName, pval = data_subset$P, id = data_subset$id),
  plink_bin = "clump/plink_win64_20231211/plink.exe",
  bfile = "clump/1kg.v3/EUR"
)
viewn(data_clumped)
# 转换并合并数据
data_clumped <- data_clumped %>%
  rename(SNP = rsid, pval.exposure = pval, id.exposure = id) %>%
  left_join(data_for_clump, by = c("SNP", "pval.exposure", "id.exposure"))

# 写入处理后的文件
output_file <- paste0("data_", file_num, "_clumped_8_F.csv")
write.csv(data_clumped_F, output_file, row.names = FALSE)
