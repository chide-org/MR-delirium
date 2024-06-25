library(TwoSampleMR)
library(LDlinkR)
library(data.table)
library(tidyverse)
homePath <- "C:/code/A-MR"
outcomeDataPath <- "C:/code/A-MR/using_outcome.csv"
tempPath <- "C:/code/A-MR/exposure/mydata/1-temp"
clumpPath <- "C:/code/A-MR/exposure/mydata/2-clump"
FPath <- "C:/code/A-MR/exposure/mydata/3-F"
harmonPath <- "C:/code/A-MR/exposure/mydata/4-harmon"
pressoPath <- "C:/code/A-MR/exposure/mydata/5-presso"
outPath <- "C:/code/A-MR/exposure/mydata/6-out"
sumPath <- "C:/code/A-MR/exposure/mydata/7-sum"

outcomedat <- fread(outcomeDataPath, header = T, nThread = 8)



# clump
setwd(tempPath)
filenames <- list.files(".")
filename <- filenames[13]
filename

{
    setwd(tempPath)
    options(ieugwasr_api = "gwas-api.mrcieu.ac.uk/")
    exposure <- read_exposure_data(
        filename,
        clump = F,
        sep = ",",
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
        samplesize_col = "N",
        min_pval = 0, # less than this value will be set to this value
        log_pval = FALSE,
    )
    exposure_clumped <- clump_data(
        exposure,
        clump_kb = 10000,
        clump_r2 = 0.001,
        clump_p1 = 1,
        clump_p2 = 1,
        pop = "EUR",
        bfile = NULL,
        plink_bin = NULL
    )
    setwd(clumpPath)
    fwrite(exposure_clumped, file = filename, row.names = F)

    # F>10
    EXP <- exposure_clumped
    if (sum(is.na(EXP$eaf.exposure)) > 0) {
        print(paste0("eaf.exposure has NA in ", filename))
        exp_f <- EXP
    } else {
        EXP$R2 <- (2 * EXP$beta.exposure * EXP$beta.exposure * EXP$eaf.exposure * (1 - EXP$eaf.exposure) / (2 * EXP$beta.exposure * EXP$beta.exposure * EXP$eaf.exposure * (1 - EXP$eaf.exposure) + 2 * EXP$se.exposure * EXP$se.exposure * EXP$samplesize.exposure * EXP$eaf.exposure * (1 - EXP$eaf.exposure))) # 计算R2
        EXP$F_value <- EXP$R2 * (EXP$samplesize.exposure - 2) / (1 - EXP$R2) # 计算F检验值
        exp_f <- EXP %>% filter(F_value > 10)
    }
    setwd(FPath)
    fwrite(exp_f, file = filename, row.names = F)

    # proxy
    missing_IVs <- exp_f$SNP[!(exp_f$SNP %in% outcomedat$SNP)]


    # Replacing missing instruments from outcome GWAS with proxies
    if (length(missing_IVs) == 0) {
        print("All exposure IVs found in outcome GWAS.")
    } else {
        print(paste0("Some exposure IVs missing from outcome GWAS: ", as.character(length(missing_IVs))))
        print("List of IVs missing from outcome GWAS:")
        for (i in 1:length(missing_IVs)) {
            print(paste0(missing_IVs[i]))
        }
        for (i in 1:length(missing_IVs)) {
            proxies <- LDproxy(snp = missing_IVs[i], pop = "EUR", r2d = "r2", token = "b6dc5de3998f", file = FALSE) # 这里的token需要在https://ldlink.nih.gov/?tab=apiaccess注册获取，这个是我的
            proxies <- proxies[proxies$R2 > 0.8, ]
            proxy_present <- FALSE

            tryCatch(
                {
                    if (length(proxies$RS_Number) == 0) {
                        print(paste0("No proxy SNP available for ", missing_IVs[i]))
                    } else {
                        for (j in 1:length(proxies$RS_Number)) {
                            proxy_present <- proxies$RS_Number[j] %in% outcomedat$SNP
                            if (proxy_present) {
                                proxy_SNP <- proxies$RS_Number[j]
                                proxy_SNP_allele_1 <- str_sub(proxies$Alleles[j], 2, 2)
                                proxy_SNP_allele_2 <- str_sub(proxies$Alleles[j], 4, 4)
                                original_SNP_allele_1 <- str_sub(proxies$Alleles[1], 2, 2)
                                original_SNP_allele_2 <- str_sub(proxies$Alleles[1], 4, 4)
                                break
                            }
                        }
                    }
                },
                error = function(e) {
                    # Handle the error
                    print(paste0("Error: ", conditionMessage(e)))
                }
            )

            if (proxy_present == TRUE) {
                print(paste0("Proxy SNP found. ", missing_IVs[i], " replaced with ", proxy_SNP))
                proxy_row <- outcomedat[outcomedat$SNP == proxy_SNP, ]
                proxy_row$SNP <- missing_IVs[i]
                if (outcomedat[outcomedat$SNP == proxy_SNP, "effect_allele.outcome"] == proxy_SNP_allele_1) proxy_row$effect_allele.outcome <- original_SNP_allele_1
                if (outcomedat[outcomedat$SNP == proxy_SNP, "effect_allele.outcome"] == proxy_SNP_allele_2) proxy_row$effect_allele.outcome <- original_SNP_allele_2
                if (outcomedat[outcomedat$SNP == proxy_SNP, "other_allele.outcome"] == proxy_SNP_allele_1) proxy_row$other_allele.outcome <- original_SNP_allele_1
                if (outcomedat[outcomedat$SNP == proxy_SNP, "other_allele.outcome"] == proxy_SNP_allele_2) proxy_row$other_allele.outcome <- original_SNP_allele_2
                outcomedat <- rbind(outcomedat, proxy_row)
            }

            if (proxy_present == FALSE) {
                print(paste0("No proxy SNP available for ", missing_IVs[i], " in outcome GWAS."))
            }
        }
    }


    # Harmon
    exp_h <- harmonise_data(
        exposure_dat = exp_f,
        outcome_dat = outcomedat,
        action = 2
    )
    setwd(harmonPath)
    fwrite(exp_h, filename, row.names = FALSE)


    # presso
    tryCatch(
        {
            presso <- run_mr_presso(exp_h)
            setwd(pressoPath)
            fwrite(
                presso[[1]]$`Main MR results`,
                paste0("Main_", filename, sep = "")
            )
            fwrite(
                presso[[1]]$`MR-PRESSO results`$`Global Test`,
                paste0("Global_", filename, sep = "")
            )
            fwrite(
                presso[[1]]$`MR-PRESSO results`$`Outlier Test`,
                paste0("Outlier_", filename, sep = "")
            )
            print(paste0("Presso for ", filename, " done!"))
        },
        error = function(e) {
            print(paste0("Error: ", filename, e))
            warning(paste0("Error: ", filename, e))
        }
    )


    # mr
    setwd(outPath)
    mr <- mr(exp_h)
    fwrite(mr, file = paste0("mr_", filename, sep = ""), row.names = F)
    tryCatch(
        {
            hetero <- mr_heterogeneity(exp_h)
            fwrite(hetero, file = paste0("hetero_", filename, sep = ""), row.names = F)
        },
        error = function(e) {
            print(paste0("Error: ", filename, e))
            warning(paste0("Error: ", filename, e))
        }
    )
    tryCatch(
        {
            pleio <- mr_pleiotropy_test(exp_h)
            fwrite(pleio, file = paste0("pleio_", filename, sep = ""), row.names = F)
        },
        error = function(e) {
            print(paste0("Error: ", filename, e))
            warning(paste0("Error: ", filename, e))
        }
    )


    # plot
    p1 <- mr_scatter_plot(
        mr_results = mr(
            exp_h,
            method_list = c(
                "mr_ivw",
                "mr_weighted_median",
                "mr_egger_regression"
            )
        ), exp_h
    )
    pdf(file = paste0("scatter_", filename, ".pdf", sep = ""), height = 8, width = 12)
    print(p1)
    dev.off()

    tryCatch(
        {
            p2 <- mr_funnel_plot(singlesnp_results = mr_singlesnp(exp_h))
            pdf(file = paste0("funnel_", filename, ".pdf", sep = ""), height = 8, width = 12)
            print(p2)
            dev.off()
        },
        error = function(e) {
            print(paste0("Error: ", filename, e))
            warning(paste0("Error: ", filename, e))
        }
    )
    tryCatch(
        {
            p3 <- mr_leaveoneout_plot(
                leaveoneout_results =
                    mr_leaveoneout(exp_h)
            )
            pdf(file = paste0("leave_", filename, ".pdf", sep = ""), height = 8, width = 12)
            print(p3)
            dev.off()
        },
        error = function(e) {
            print(paste0("Error: ", filename, e))
            warning(paste0("Error: ", filename, e))
        }
    )
}
