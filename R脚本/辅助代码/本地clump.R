# 本地去除连锁不平衡的SNP
library(plinkbinr)
library(ieugwasr)
library(data.table)
clumptoolPath <- "D:/A-MR/clump_tool"
tempPath <- "D:/A-MR/exposure/mydata/temp"
clumpPath <- "D:/A-MR/exposure/mydata/clump"

setwd(tempPath)
list <- list.files(".")
list
for (file in list) {
  tryCatch(
    {
      setwd(tempPath)
      dat <- fread(file, header = T)
      setwd(clumptoolPath)
      getwd()
      # 教程：https://mp.weixin.qq.com/s/TDEZYlp2AHHA8Yqld8GKHA
      # 输入待clump的数据
      dat <- dat %>%
        rename(rsid = SNP, pval = p) %>%
        as_tibble()
      data_clumped <- ld_clump(
        dat,
        clump_kb = 10000,
        clump_r2 = 0.001,
        clump_p = 0.99,
        pop = "EUR",
        plink_bin = "plink_win64_20231211/plink.exe",
        bfile = "1kg.v3/EUR"
      )
      data_clumped <- data_clumped %>% rename(SNP = rsid, p = pval)
      # 输出clump结果
      setwd(clumpPath)
      fwrite(data_clumped, file = paste0("clump_", file, sep = ""))
    },
    error = function(e) {
      print(paste0("Error: ", file, "\n", e))
      return()
    }
  )
}
