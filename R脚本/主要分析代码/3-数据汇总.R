library(tidyverse)
library(data.table)
library(TwoSampleMR)
homePath <- "C:/code/A-MR"
outcomeDataPath <- "C:/code/A-MR/using_outcome.csv"
tempPath <- "C:/code/A-MR/exposure/mydata/1-temp"
clumpPath <- "C:/code/A-MR/exposure/mydata/2-clump"
FPath <- "C:/code/A-MR/exposure/mydata/3-F"
harmonPath <- "C:/code/A-MR/exposure/mydata/4-harmon"
pressoPath <- "C:/code/A-MR/exposure/mydata/5-presso"
outPath <- "C:/code/A-MR/exposure/mydata/6-out"
sumPath <- "C:/code/A-MR/exposure/mydata/7-sum"
# 导入表型名单
setwd(tempPath)
filenames <- list.files(".", full.names = F)
filenames


# 创建一个空表格，并指定列名
mr_output <- data.frame(
  Exposure = character(),
  ID.exposure = character(),
  N.samples = character(),
  nsnp = numeric()
)
# 其他列根据method列自动生成
# cycle for mr_output
for (filename in filenames) {
  tryCatch(
    {
      # 读取表型数据
      setwd(outPath)
      mrdat <- fread(paste0("mr_", filename), header = T)
      mrdat <- generate_odds_ratios(mrdat)
      Exposure <- mrdat$exposure[1]
      ID.exposure <- mrdat$id.exposure[1]
      setwd(tempPath)
      N <- fread(filename, header = T)
      if (length(unique(N$N)) == 1) {
        N.samples <- N$N[1]
      } else {
        N.samples <- "not one!"
      }
      nsnp <- mrdat$nsnp[1]


      # 生成新列名
      methods <- unique(mrdat$method)
      new_columns <- c()
      for (method in methods) {
        new_columns <- c(
          new_columns,
          paste0(
            c(
              "pval",
              "lo_ci",
              "up_ci",
              "b",
              "se",
              "or",
              "or_lci95",
              "or_uci95"
            ),
            "_",
            method
          )
        )
      }
      # 将方法名称转换为列名
      temp_data <- mrdat %>%
        pivot_wider(
          names_from = method,
          values_from = c(
            "pval",
            "lo_ci",
            "up_ci",
            "b",
            "se",
            "or",
            "or_lci95",
            "or_uci95"
          )
        ) %>%
        select(-outcome, -exposure, -id.outcome)

      # 定义一个函数，用于替换每列的NA值为该列的唯一非NA值，用于扩展合并不同行为一行
      replace_na_with_unique <- function(x) {
        unique_values <- na.omit(x) # 获取非NA值
        if (length(unique_values) == 0) {
          # 如果没有非NA值
          return(NA)
        } else {
          return(replace(x, is.na(x), unique_values[1])) # 将NA值替换为唯一非NA值
        }
      }

      # 压缩为一行
      temp_data <- temp_data %>%
        mutate(across(everything(), replace_na_with_unique)) %>%
        summarise_all(list(~ first(.)))
      temp_data <- temp_data %>%
        mutate(
          Exposure = Exposure,
          ID.exposure = ID.exposure,
          N.samples = N.samples,
          nsnp = nsnp
        ) %>%
        dplyr::select("Exposure", "ID.exposure", "N.samples", "nsnp", all_of(new_columns))

      # 合并结果
      # 获取所有可能的列名
      all_columns <- union(names(mr_output), names(temp_data))
      if (nrow(mr_output)) {
        mr_output[, setdiff(all_columns, names(mr_output))] <- NA
        mr_output <- mr_output[, all_columns]
      }
      temp_data[, setdiff(all_columns, names(temp_data))] <- NA
      temp_data <- temp_data[, all_columns]
      mr_output <- rbind(mr_output, temp_data)
    },
    error = function(e) {
      print(paste("Error: ", filename, e))
    }
  )
}
View(mr_output)
setwd(sumPath)
write.csv(mr_output, "mr_sum.csv", row.names = F)


# 创建一个空表格，并指定列名
ph_output <- data.frame(
  Exposure = character(),
  ID.exposure = character(),
  N.samples = character(),
  nsnp = numeric(),
  egger_intercept = numeric(),
  se_intercept = numeric(),
  pval_intercept = numeric(),
  Q_pval_ivw = numeric(),
  Q_pval_egger = numeric()
)
# cycle for ph_output

for (filename in filenames) {
  tryCatch(
    {
      setwd(outPath)
      pleio <- fread(paste0("pleio_", filename), header = T)
      hetero <-
        fread(paste0("hetero_", filename), header = T)

      ivw <-
        filter(hetero, method == "Inverse variance weighted")$Q_pval
      egger <- filter(hetero, method == "MR Egger")$Q_pval
      setwd(tempPath)
      N <- fread(filename, header = T)
      if (length(unique(N$N)) == 1) {
        N.samples <- N$N[1]
      } else {
        N.samples <- "not one!"
      }
      temp_data <- pleio %>%
        select(
          Exposure = exposure,
          ID.exposure = id.exposure,
          egger_intercept = egger_intercept,
          se_intercept = se,
          pval_intercept = pval
        ) %>%
        mutate(
          N.samples = N.samples,
          Q_pval_ivw = ivw,
          Q_pval_egger = egger
        )
      ph_output <- rbind(ph_output, temp_data)
    },
    error = function(e) {
      print(paste("Error: ", filename, e))
    }
  )
}
View(ph_output)
setwd(sumPath)
write.csv(ph_output, "ph_sum.csv", row.names = F)


# 创建一个空表格，并指定列名
confounder_output <- data.frame(
  Exposure = character(),
  ID.exposure = character(),
  SNP = character(),
  Confounder = character(),
  PMID = character(),
  P.value = numeric(),
  Detection.method = character()
)
# cycle for confounder_output
setwd("D:/A-MR/exposure")
getwd()
confuse <- fread("Confounder_SNP.csv", header = T) %>% select(-X)

setwd("D:/A-MR/exposure/clump")
for (filename in filenames) {
  tryCatch(
    {
      snp <-
        fread(paste0(filename), header = T) %>% select("SNP")
      setwd(manyDataPath)
      mrdat <- fread(paste0("mr_", filename), header = T)
      setwd("D:/A-MR/exposure/clump")
      Exposure <- mrdat$exposure[1]
      ID.exposure <- mrdat$id.exposure[1]
      nsnp <- mrdat$nsnp[1]
      temp_data <-
        confuse[confuse$SNP %in% snp$SNP, ] %>% mutate(Exposure, ID.exposure, nsnp, .before = 1)
      confounder_output <- rbind(confounder_output, temp_data)
    },
    error = function(e) {
      print(paste0("Error in ", filename, e))
    }
  )
}
View(confounder_output)
setwd("D:/A-MR/exposure")
write.csv(confounder_output, "confounder_output.csv", row.names = F)
