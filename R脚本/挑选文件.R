library(tidyverse)

setwd("D:/A-MR/exposure/confused")
list <- fread("../有去除混淆snp的表型.csv", header = F)
allDir <- "../confused/"
files <- list.files(allDir, full.names = F)
files

# data deal
list <- list$V1 %>% substring(8)
list

for (file in files) {
    tryCatch(
        {
            if (substring(file, 9) %in% list) {
                file.rename(file, paste0("../需要重分析的数据/", file, sep = ""))
            }
        },
        error = function(e) {
            print(e)
        }
    )
}
