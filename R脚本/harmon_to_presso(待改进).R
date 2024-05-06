library(tidyverse)
library(stringr)
library(data.table)

setwd("D:/A-MR")
harmonDir <- "exposure/harmon/"
harmonFiles <- list.files(harmonDir)

setwd("D:/A-MR/exposure/harmon")
harmonFiles <- c("harmon_TSH.csv", "harmon_FT4.csv")


for (harmonFile in harmonFiles) {
    tryCatch(
        {
            dat_harmon <- fread(harmonFile)
            dat_harmon <- as.data.frame(dat_harmon)
            presso <- run_mr_presso(dat_harmon, NbDistribution = 5000)

            dat_presso <- dat_harmon[!(dat_harmon$SNP %in% presso$Outliers)] # 需要测试

            write.csv(dat_presso, file = paste0("out/", "presso_", substring(harmonFile, 8), sep = ""))
            print(paste0("Presso for ", harmonFile, " done!"))
        },
        error = function(e) {
            print(paste0("Error: ", conditionMessage(e)))
        }
    )
}

typeof(dat_harmon)
as.data.frame(dat_harmon)

for (harmonFile in harmonFiles) {
    dat_harmon <- fread(paste0(harmonDir, harmonFile, sep = ""))
    dat_harmon$outcome <- "delirium"
    write.csv(dat_harmon, file = paste0(harmonDir, harmonFile, sep = ""))
}
