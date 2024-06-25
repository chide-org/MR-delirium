library(tidyverse)
library(stringr)
library(data.table)

needPressos <- none

# main
for (needPresso in needPressos) {
  tryCatch({
    setwd("D:/A-MR/exposure/harmon")
    dat_harmon <-
      read.csv(paste0("harmon_", needPresso, ".csv", sep = ""), header = T)
    #dat_harmon <- as.data.frame(dat_harmon)
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
    #dat_presso <- dat_harmon[!(dat_harmon$SNP %in% presso$Outliers)] # 需要测试
    
    #write.csv(dat_presso, file = paste0("out/", "presso_", substring(needPresso, 8), sep = ""))
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
