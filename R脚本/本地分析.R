cycle(){
    outcomePath <- "outcome/delirium.csv"
    exposure_dirPath <- "exposure"
    # 处理结局
    outcome <- fread(outcomePath)
    viewn(outcome)
    nrow(outcome)
    excel_open(outcome)
    outcome <- excel_close(outcome)
    # 批量处理暴露
    names <- list.files(exposure_dirPath)
    paths <- paste0(exposure_dirPath, "/", names, sep = "")
    paths
    i <- 1
    # 空白行，初始化all_mr
    empty_row <- data.frame(matrix(NA, nrow = 1, ncol = 9))
    colnames(empty_row) <- c("id.exposure","id.outcome","outcome","exposure","method","nsnp","b","se","pval")
    all_mr <- empty_row
    # 输出
    write.csv(all_mr, file = "all_mr.csv")
    # 循环
    for(path in paths){
        # exposure
        exposure <- fread(path)
        dat <- merge(exposure, outcome, by.x = "SNP", by.y = "SNP")
        write.csv(dat, file = "temp.csv")
        dat_read <- read_outcome_data(
        snps = exposure$SNP,
        filename = "temp.csv",
        sep = ",",
        snp_col = "SNP",
        pval_col = "Pvalue",
        effect_allele_col = "Effect_allele",
        other_allele_col = "Non_Effect_allele",
        beta_col = "Beta",
        se_col = "SE",
        eaf_col = "EAF",
        )
        file.remove("temp.csv")
        # harmon
        dat_harmon <- harmonise_data(
            exposure_dat = exposure,
            outcome_dat = dat_read
        )
        write.csv(dat_harmon, file = paste0(path,"_harmon.csv"))
        # mr
        mr <- mr(dat_harmon)
        all_mr <- rbind(all_mr, empty_row, mr)
        # 代码还需要添加代理SNP的处理和敏感性分析
    }
}

single(){
    outcomePath <- "outcome/delirium.csv"
    exposurePath <- "exposure/（多）ukb-a-158.csv"
    exposure_name <- "测试"
    # 处理结局(批量处理执行一次即可)
    outcome <- fread(outcomePath)
    colnames(outcome)
    viewn(outcome)
    nrow(outcome)
    excel_open(outcome)
    #！！！修改列名并保存后再运行下一条
    outcome <- excel_close(outcome)
    # 处理暴露
    exposure <- fread(exposurePath)
    dat <- merge(exposure,outcome,by.x = "SNP",by.y="SNP")
    write.csv(dat, file = "temp.csv")
    dat_read <- read_outcome_data(
        snps = exposure$SNP,
        filename = "temp.csv",
        sep = ",",
        snp_col = "SNP",
        pval_col = "Pvalue",
        effect_allele_col = "Effect_allele",
        other_allele_col = "Non_Effect_allele",
        beta_col = "Beta",
        se_col = "SE",
        eaf_col = "EAF",
        id_col = "id",
        phenotype_col = "phenotype"
        )
    file.remove("temp.csv")
    # harmon
    dat_harmon <- harmonise_data(
        exposure_dat = exposure,
        outcome_dat = dat_read
    )
    write.csv(dat_harmon,file =paste0(exposure_name,"_harmon.csv",sep = ""))
    # mr
    mr <- mr(dat_harmon)
    mr$outcome <- "delirium"
    mr$id.outcome <- "delirium"
    mr$exposure <- exposure_name
    mr$id.exposure <- exposure_name
    View(mr)
    write.csv(mr,file = paste0(exposure_name,"_mr.csv",sep = ""))
}

#敏感性分析(须在R中输出图片)
mr_scatter_plot(
  mr_results = mr(
    dat_harmon,
    method_list = c("mr_ivw",
                  "mr_weighted_median",
                  "mr_egger_regression"
                             )
    ),
  dat_harmon
)

dat_harmon_hetero <- mr_heterogeneity(dat_harmon)
View(dat_harmon_hetero)
write.csv(dat_harmon_hetero,file = paste0(exposure_name,"_hetero.csv",sep = ""))




mr_funnel_plot(singlesnp_results = mr_singlesnp(dat_harmon))

pleio <- mr_pleiotropy_test(dat_harmon)
View(pleio)
write.csv(pleio,file = paste0(exposure_name,"_pleio.csv",sep = ""))

mr_leaveoneout_plot(
  leaveoneout_results =
    mr_leaveoneout(dat_harmon)
)
