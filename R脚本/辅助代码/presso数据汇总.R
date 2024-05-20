library(tidyverse)

# 导入表型名单
dir <- "D:/A-MR/exposure/clump"
phenotype_list <- list.files(dir, full.names = F)
phenotype_list <- substr(phenotype_list, 1, nchar(phenotype_list) - 4)
phenotype_list

# 定义一个函数，用于替换每列的NA值为该列的唯一非NA值
replace_na_with_unique <- function(x) {
  unique_values <- na.omit(x) # 获取非NA值
  if (length(unique_values) == 0) {
    # 如果没有非NA值
    return(NA)
  } else {
    return(replace(x, is.na(x), unique_values[1])) # 将NA值替换为唯一非NA值
  }
}

setwd("D:/A-MR/exposure/无混淆presso结果")
presso_output <- data.frame()
none <- c()
# cycle
for (phenotype in phenotype_list) {
  tryCatch({
    # 读取表型数据
    global <-
      read.csv(paste0("mr_presso_Globalharmon_", phenotype,".csv", sep = ""), header = T)

    temp_data <-
      read.csv(paste0("mr_presso_Mainharmon_", phenotype,".csv", sep = ""),
               header = T) %>%
      pivot_wider(
        names_from = MR.Analysis,
        values_from = c("Causal.Estimate", "Sd", "T.stat", "P.value")
      )%>%
      mutate(RSSobs = global$RSSobs[1],Pvalue = global$Pvalue[1])%>%
      mutate(across(everything(), replace_na_with_unique)) %>%
      summarise_all(list( ~ first(.))) %>%
      mutate(ID.exposure = phenotype,.before = 1) %>%
      select(-X,-Exposure)


    presso_output <- rbind(presso_output,temp_data)
  },
  error = function(e) {
    print(paste("Error: ", phenotype))
    none <<- c(none,phenotype)
  })
}
View(presso_output)
write.csv(presso_output, "无混淆_presso_output.csv", row.names = F)
