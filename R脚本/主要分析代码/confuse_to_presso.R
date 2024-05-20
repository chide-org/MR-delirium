library(tidyverse)
library(stringr)
# 导入表型名单
dir <- "D:/A-MR/exposure/clump"
needPressos <- list.files(dir, full.names = F)
needPressos <- substr(needPressos, 1, nchar(needPressos) - 4)
needPressos

# main
for (needPresso in needPressos) {
  tryCatch({
    setwd("D:/A-MR/exposure/confused")
    dat_harmon <-
      read.csv(paste0("confuse_", needPresso, ".csv", sep = ""), header = T)
    presso <- run_mr_presso(dat_harmon)
    setwd("D:/A-MR/exposure/presso")
    write.csv(
      presso[[1]]$`Main MR results`,
      paste0("mr_presso_Mainharmon_", needPresso, ".csv",sep = "")
    )
    write.csv(
      presso[[1]]$`MR-PRESSO results`$`Global Test`,
      paste0("mr_presso_Globalharmon_", needPresso, ".csv",sep = "")
    )
    write.csv(
      presso[[1]]$`MR-PRESSO results`$`Outlier Test`,
      paste0("mr_presso_Outlierharmon_", needPresso, ".csv",sep = "")
    )
    print(paste0("Presso for ", needPresso, " done!"))
  },
  error = function(e) {
    print(paste0("Error: ",needPresso))
  })
}

# co
for (needPresso in needPressos) {
  dat_harmon <- fread(paste0(harmonDir, needPresso, sep = ""))
  dat_harmon$outcome <- "delirium"
  write.csv(dat_harmon, file = paste0(harmonDir, needPresso, sep = ""))
}
