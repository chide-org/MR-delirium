library(tidyverse)
setwd("D:/A-MR/exposure/需要手动汇总的暴露/有混淆")

# 导入表型名单
dir <- "D:/A-MR/exposure/clump"
phenotype_list <- list.files(".", full.names = F)
phenotype_list <- phenotype_list[grep("mr_*",phenotype_list)]
phenotype_list

deff_mr <- data.frame()

# cycle
for (phenotype in phenotype_list) {
  tryCatch({
    # 读取表型数据
    mrdat <-read.csv(phenotype, header = T)
    mrdat <- generate_odds_ratios(mrdat)
    deff_mr <-rbind(deff_mr,mrdat)
  },
  error = function(e) {
    print(paste("Error: ", phenotype))
  })
}
View(deff_mr)
write.csv(deff_mr, "all_deff_mr.csv", row.names = F)


