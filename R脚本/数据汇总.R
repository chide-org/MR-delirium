library(tidyverse)
setwd("D:/A-MR")
# 创建一个空表格，并指定列名
mr_output <- data.frame(Exposure = character(),
                        ID.exposure = character(),
                        nsnp = numeric())
# 其他列根据method列自动生成

# 创建一个空表格，并指定列名
ph_output <- data.frame(
  Exposure = character(),
  ID.exposure = character(),
  nsnp = numeric(),
  egger_intercept = numeric(),
  se_intercept = numeric(),
  pval_intercept = numeric(),
  Q_pval_ivw = numeric(),
  Q_pval_egger = numeric()
)

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

# 导入表型名单
dir <- "D:/A-MR/exposure/clump"
phenotype_list <- list.files(dir, full.names = F)
phenotype_list <-
  substr(phenotype_list, 1, nchar(phenotype_list) - 4)
phenotype_list
#manyDataPath <-"D:/A-MR/exposure/含混淆因素-csv"
manyDataPath <-"D:/A-MR/exposure/无混淆-csv"



# cycle for mr_output
setwd(manyDataPath)
for (phenotype in phenotype_list) {
  tryCatch({
    # 读取表型数据
    mrdat <-read.csv(paste0("mr_", phenotype, ".csv"), header = T)
 #   setwd("D:/A-MR/exposure/需要手动汇总的暴露")
  #  mrdat <- read.csv("mr_（多）ukb-b-13532.csv",header = T)
    mrdat <- generate_odds_ratios(mrdat)
    Exposure <- mrdat$exposure[1]
    ID.exposure <- mrdat$id.exposure[1]
    nsnp <- mrdat$nsnp[1]
    
    # 生成新列名
    methods <- unique(mrdat$method)
    new_columns <- c()
    for (method in methods) {
      new_columns <- c(new_columns,
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
                       ))
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
    
    # 压缩为一行
    temp_data <- temp_data %>%
      mutate(across(everything(), replace_na_with_unique)) %>%
      summarise_all(list(~ first(.)))
    temp_data <- temp_data %>%
      mutate(Exposure = Exposure,
             ID.exposure = ID.exposure,
             nsnp = nsnp) %>%
      dplyr::select("Exposure", "ID.exposure", "nsnp", all_of(new_columns))
    
    # 合并结果
    
    mr_output <- rbind(mr_output, temp_data)
  },
  error = function(e) {
    print(paste("Error: ", phenotype))
  })
}
View(mr_output)
setwd("D:/A-MR/exposure")
write.csv(mr_output, "无混淆_mr_output.csv", row.names = F)



# cycle for ph_output
setwd(manyDataPath)
for (phenotype in phenotype_list) {
  tryCatch({
    pleio <- read.csv(paste0("pleio_", phenotype, ".csv"), header = T)
    hetero <-
      read.csv(paste0("hetero_", phenotype, ".csv"), header = T)
    
    ivw <-
      filter(hetero, method == "Inverse variance weighted")$Q_pval
    egger <- filter(hetero, method == "MR Egger")$Q_pval
    
    temp_data <- pleio %>%
      select(
        Exposure = exposure,
        ID.exposure = id.exposure,
        egger_intercept,
        intercept_se = se,
        intercept_pval = pval
      ) %>%
      mutate(Q_pval_ivw = ivw,
             Q_pval_egger = egger)
    ph_output <- rbind(ph_output, temp_data)
  },
  error = function(e) {
    print(paste("Error: ", phenotype))
  })
}
View(ph_output)
setwd("D:/A-MR/exposure")
write.csv(ph_output, "无混淆_ph_output.csv", row.names = F)



# cycle for confounder_output
setwd("D:/A-MR/exposure")
getwd()
confuse <- read.csv("Confounder_SNP.csv", header = T) %>% select(-X)

setwd("D:/A-MR/exposure/clump")
for (phenotype in phenotype_list) {
  tryCatch({
    snp <-
      read.csv(paste0(phenotype, ".csv"), header = T) %>% select(SNP)
    temp_data <-
      confuse[confuse$SNP %in% snp$SNP, ] %>% mutate(Exposure = phenotype, .before = 1)
    confounder_output <- rbind(confounder_output, temp_data)
  },
  error = function(e) {
    print(paste0("Error in ", phenotype, e))
  })
}
View(confounder_output)
setwd("D:/A-MR/exposure")
write.csv(confounder_output, "confounder_output.csv", row.names = F)
