library(TwoSampleMR)
library(data.table)
library(tidyverse)
# 设置工作路径
exposurePath <- "D:/A-MR/exposure/mydata" # '\'要手动换成'/'
TwoSampleMRPath <- "D:/A-MR/library/TwoSampleMR" # '\'要手动换成'/'
tempPath <- "D:/A-MR/exposure/mydata/temp"

# START
setwd(exposurePath)
getwd()
# 读取原始数据,压缩包解压后读取更快,bandzip
bigdat <- fread("bip2019.gz", header = T) # !!!!!!!!
name <- "bip2019" # !!!!!!!

# 修改列名
# 必须列：SNP,EA,OA,se,p,beta
# 必要列：eaf,CHR,BP,id,phe,nca,ncon,N
excel_open(bigdat)
bigdat <- excel_close(bigdat)

# P值筛选
dat <- subset(bigdat, p < 5e-6)
nrow(dat)
# !!!检查是否有数据再运行下面的代码
rm(bigdat)


# 检查SNP列是否是否可用
temp <- grepl("^rs*", dat$SNP) # 正则表达式匹配SNP列中以rs开头的行

sum(!temp) # 不符合要求的总行数，大于0 => JUMP TO chr_to_rs.R


# 必须列：SNP,EA,OA,se,p,beta
# 必要列：eaf,CHR,BP,id,phe,nca,ncon,N
# 没有则以NA代替
colnames(dat)
dat$or <- as.numeric(dat$or)

dat <- dat %>% mutate(
  # beta = log(or),
  eaf = (f1 + f2) / 2,
  # eaf = NA,
  id = name,
  phe = name,
  # nca = NA,
  # ncon = NA,
  # N = NA,
  N = nca + ncon,
)

# 挑选所有需要的列
dat <- select(
  dat,
  c(
    "SNP", "EA", "OA", "p", "se", "beta",
    "eaf", "CHR", "BP", "id", "phe", "nca", "ncon", "N"
  )
)
viewn(dat)

# 暂存规范化的数据
setwd(tempPath)
getwd()
fwrite(dat, file = paste0(name, ".csv", sep = ""), nThread = 8)

# p < 5e-6
fwrite(dat, file = paste0(name, "_p5e-6.csv", sep = ""), nThread = 8)
