# 将待clump的数据写成文件移动到TwoSampleMR包下
setwd(TwoSampleMRPath)
getwd()
write.csv(dat, file = "temp.csv", row.names = F)

road <- system.file("temp.csv", package = "TwoSampleMR") # 获取文件路径


# 进行clump,失败则进行本地clump
dat_read <- read_exposure_data(
  filename = road,
  clump = F,
  sep = ",", #  csv文件以','作为分隔符
  
  snp_col = "SNP",
  effect_allele_col = "EA",
  other_allele_col = "OA",
  pval_col = "p",
  se_col = "se",
  beta_col = "beta",
  eaf_col = "eaf",
  chr_col = "CHR",
  pos_col = "BP",
  id_col = "id",
  phenotype_col = "phe",
  ncase_col = "nca",
  ncontrol_col = "ncon",
  samplesize_col = "N"
)

# clump失败则执行以下代码
{
  # clump失败 JUMP TO 本地clump得到dat_clumped
  dat_read <- dat_read[dat_read$SNP %in% data_clumped$rsid, ]
}

viewn(dat_read)
setwd(work)
write.csv(dat_read, file = "clump/adhd2019.csv", row.names = F)
