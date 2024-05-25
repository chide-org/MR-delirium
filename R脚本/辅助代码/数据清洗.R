library(TwoSampleMR)
library(data.table)
library(tidyverse)


  # 设置工作路径
work <- "D:/A-MR/exposure/rough_data"   # '\'要手动换成'/'
TwoSampleMR <- "D:/A-MR/library/TwoSampleMR"  # '\'要手动换成'/'
setwd(work)
getwd()


  # 读取原始数据,压缩包解压后读取更快,bandzip
dat <- fread("an2017-pgc.ed.freeze1.summarystatistics.July2017.txt.gz")


  # 检查是否有必须列SNP,A1,A2,se,p,beta(or)
  # 检查是否有可用列eaf(maf),ncase,ncontrol,samplesize,phenotype,id,CHR,BP
  # phenotype,id如果没有就自己创建，统一用表格"publication"下的暴露名
viewn(dat)


  # 检查SNP列是否是否可用
temp <- grepl("^rs*",dat$SNP)  # 正则表达式匹配SNP列中以rs开头的行
nrow(dat)  # 总行数
sum(temp)  # 符合要求的总行数
wrongdat <- dat[!temp]  # 将错误的行挑出来
View(wrongdat)  # !!!如果大部分是chr:pos号,那么需要进行chr转rs号,不能直接丢掉


  # JUMP TO chr_to_rs.R得到rightdat


  #合并
dat <- dat[temp]  # 将可用的行挑出来
nrow(dat)
nrow(rightdat)
dat <-rbind(dat,rightdat)
nrow(dat)
viewn(dat)


  # 修改列名准备补充列:SNP,effect_allele,other_allele,se,p,beta(or)
  # eaf(maf),ncase,ncontrol,samplesize,phenotype,id,CHR,BP
excel_open(dat)
dat <- excel_close(dat)
viewn(dat)



  # P值筛选
dat <- subset(dat,p<5e-8)
nrow(dat)
View(dat)  # 检查p值是不是都是小于5e-8


  # 根据beta=log(or),samplesize=ncase+ncontrol,se+p=>eaf,计算F值等添加列
dat$se <- as.numeric(dat$se)
dat <- dat %>% mutate(beta = log(se),id = "an2017",phenotype = "an2017")
viewn(dat)


dat <- select(
  dat,
  c(
    "SNP",
    "effect_allele",
    "other_allele",
    "beta",
    "se",
    "p",
    
    "CHR",
    "BP",
    "id",
    "phenotype"
  )
)
viewn(dat)



  # 将待clump的数据写成文件移动到TwoSampleMR包下
setwd(TwoSampleMR)
getwd()
write.csv(dat,file = "temp.csv",row.names = F)

road <- system.file("temp.csv",package = "TwoSampleMR")   # 获取文件路径


  # 进行clump,失败则进行本地clump
dat <- read_exposure_data(
  
  filename = road,
  clump = F,
  sep = ",",  #  csv文件以','作为分隔符
  
  snp_col = "SNP",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  beta_col = "beta",
  se_col = "se",
  pval_col = "p",
  
  chr_col = "CHR",
  pos_col = "BP",
  phenotype_col = "phenotype",
  id_col = "id"
)
viewn(dat)


  #clump失败 JUMP TO 本地clump得到dat_clumped


dat <- dat[dat$rsid%in%data_clumped$rsid]
sum(!dat$rsid%in%data_clumped$rsid)
length(unique(data_clumped$rsid))
nrow(dat_read)
viewn(dat_read)
setwd(work)

write.csv(dat,file = "clump/an2017.csv")





















