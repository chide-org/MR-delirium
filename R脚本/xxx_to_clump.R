# 读取暴露数据
exposure <- fread("exposure/T2D.csv")
viewn(exposure)
excel_open(exposure)
exposure <- excel_close(exposure)
colnames(exposure)
nrow(exposure)

# p
exposure <- subset(exposure, Pvalue < 5e-8)
nrow(exposure)
write.csv(exposure, file = "Ptemp.csv")

# 将bmi文件移动到TowSampleMR包下
success <-
    file.rename(
        "./Ptemp.csv",
        "./library/TwoSampleMR/Ptemp.csv"
    )
if (success) {
    print("文件移动成功！")
} else {
    print("文件未移动成功")
}
road <- system.file("Ptemp.csv", package = "TwoSampleMR")

# clump

exposure_clumped <- read_exposure_data(
    filename = road,
    sep = ",",
    snp_col = "SNP",
    pval_col = "Pvalue",
    effect_allele_col = "Effect_allele",
    other_allele_col = "Non_Effect_allele",
    beta_col = "Beta",
    se_col = "SE",
    eaf_col = "EAF",
    clump = TRUE
)
nrow(exposure_clumped)
write.csv(exposure_clumped, file = "exposure/T2D_clumped.csv")

# 在线数据循环获取
ieus <- c("ukb-a-190", "ukb-a-523", "ukb-a-76", "ukb-a-77", "ukb-b-16956", "ukb-b-17918", "ukb-b-19732", "ukb-b-20289", "ukb-b-4226", "ukb-b-9971", "ebi-a-GCST90013893", "ebi-a-GCST90013933", "ebi-a-GCST90013943", "ebi-a-GCST90018853", "ebi-a-GCST90018855", "ebi-a-GCST90018860", "ebi-a-GCST90018862", "ebi-a-GCST90018929", "ebi-a-GCST90018990", "ebi-a-GCST90029022")

for (ieu_id in ieus) {
    exposure <- extract_instruments(outcomes = ieu_id)
    write.csv(exposure, file = )
}

Sys.getenv("R_ENVIRON_USER")
user()
ieugwasr::get_opengwas_jwt()
exposure <- extract_instruments(outcomes = "ukb-a-190", access_token = ieugwasr::get_opengwas_jwt())

remotes::install_github("MRCIEU/TwoSampleMR")
remotes::install_github("MRCIEU/ieugwasr")
# 老版本
remotes::install_github("MRCIEU/TwoSampleMR@0.4.26")
remotes::install_github("MRCIEU/ieugwasr@v0.2.1")
library(TwoSampleMR)
library(ieugwasr)
exposure <- extract_instruments(outcomes = "ukb-a-190")

# 在线数据单个获取
library(TwoSampleMR)
exposure <- extract_instruments(outcomes = "prot-a-530")
write.csv(exposure, file = )
