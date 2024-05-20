library(TwoSampleMR)
library(tidyverse)
library(dplyr)
library(stringr)

setwd("D:/A-MR/")
confusePath <- "exposure/Confounder_SNP.csv"
exposure_dirPath <- "exposure/harmon"

confuse <- fread(confusePath) %>% select(SNP)
viewn(confuse)

# 批量处理暴露
names <- list.files(exposure_dirPath)
names

setwd("D:/A-MR/exposure/harmon")
getwd()

# 循环处理文件
for (file_name in names) {
    tryCatch(
        {
            presso_dat <- fread(file_name)
            # presso_dat <- presso_dat %>% select(coln)
            # write.csv(presso_dat, file = file_name)
            x <- sum(presso_dat$SNP %in% confuse$SNP)
            if (x > 0) {
                print(paste0("从", file_name, "中删除与混杂SNP", x, "个", sep = ""))
            }
            confuse_dat <- presso_dat %>% filter(!(presso_dat$SNP %in% confuse$SNP))
            write.csv(confuse_dat, file = paste0("../confused/confuse_", substring(file_name, 8), sep = ""))
        },
        error = function(e) {
            print(paste0("Error in ", file_name, e))
        }
    )
}
View(presso_dat)
coln <- colnames(presso_dat)
