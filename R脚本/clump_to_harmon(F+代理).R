library(TwoSampleMR)
library(LDlinkR)
library(tidyverse)
library(dplyr)
library(stringr)
# file_names_subset <- file_names
setwd("D:/A-MR/")
outcomePath <- "outcome/no_empty_delirium.csv"
exposure_dirPath <- "exposure"
# 处理结局
outcome <- fread(outcomePath)
viewn(outcome)
excel_open(outcome)
outcome <- excel_close(outcome)
# 批量处理暴露
names <- list.files(exposure_dirPath)
names <- c("TSH.tsv", "FT4.tsv")
setwd("D:/A-MR/exposure/clump")
getwd()

# fwrite(outcome, "D:/A-MR/outcome/outcome_delirium.csv")
# file_name <- "FT4.tsv"
# 循环处理文件
out_dat <- outcome
# Loop through each file
for (file_name in names) {
    tryCatch(
        {
            # Read the CSV file
            EXP <- read.delim(file_name)

            if (sum(is.na(EXP$eaf.exposure)) > 0) {
                print(paste0("eaf.exposure has NA in ", file_name))
                exp_dat_clumped <- EXP
            } else {
                # 计算F检验值
                EXP$R2 <- (2 * EXP$beta.exposure * EXP$beta.exposure * EXP$eaf.exposure * (1 - EXP$eaf.exposure) / (2 * EXP$beta.exposure * EXP$beta.exposure * EXP$eaf.exposure * (1 - EXP$eaf.exposure) + 2 * EXP$se.exposure * EXP$se.exposure * EXP$samplesize.exposure * EXP$eaf.exposure * (1 - EXP$eaf.exposure))) # 计算R2
                EXP$F_value <- EXP$R2 * (EXP$samplesize.exposure - 2) / (1 - EXP$R2) # 计算F检验值

                exp_dat_clumped <- EXP %>% filter(F_value > 10)
            }

            # Identifying & printing exposure instruments missing from outcome GWAS
            missing_IVs <- exp_dat_clumped$SNP[!(exp_dat_clumped$SNP %in% out_dat$SNP)]
            print(paste0("Number of IVs missing from outcome GWAS: ", as.character(length(missing_IVs))))
            print("List of IVs missing from outcome GWAS:")
            for (i in 1:length(missing_IVs)) {
                print(paste0(missing_IVs[i]))
            }

            # Replacing missing instruments from outcome GWAS with proxies
            if (length(missing_IVs) == 0) {
                print("All exposure IVs found in outcome GWAS.")
            } else {
                print("Some exposure IVs missing from outcome GWAS.")
                out_full <- outcome

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
                                    proxy_present <- proxies$RS_Number[j] %in% out_full$SNP

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
                        proxy_row <- out_dat[1, ]
                        proxy_row$SNP <- missing_IVs[i]
                        proxy_row$beta.outcome <- as.numeric(out_full[out_full$SNP == proxy_SNP, "beta.outcome"])
                        proxy_row$se.outcome <- as.numeric(out_full[out_full$SNP == proxy_SNP, "se.outcome"])
                        if (out_full[out_full$SNP == proxy_SNP, "effect_allele.outcome"] == proxy_SNP_allele_1) proxy_row$effect_allele.outcome <- original_SNP_allele_1
                        if (out_full[out_full$SNP == proxy_SNP, "effect_allele.outcome"] == proxy_SNP_allele_2) proxy_row$effect_allele.outcome <- original_SNP_allele_2
                        if (out_full[out_full$SNP == proxy_SNP, "other_allele.outcome"] == proxy_SNP_allele_1) proxy_row$other_allele.outcome <- original_SNP_allele_1
                        if (out_full[out_full$SNP == proxy_SNP, "other_allele.outcome"] == proxy_SNP_allele_2) proxy_row$other_allele.outcome <- original_SNP_allele_2
                        proxy_row$pval.outcome <- as.numeric(out_full[out_full$SNP == proxy_SNP, "pval.outcome"])
                        if ("eaf.outcome" %in% colnames(out_full)) proxy_row$eaf.outcome <- as.numeric(out_full[out_full$SNP == proxy_SNP, "eaf.outcome"])
                        out_dat <- rbind(out_dat, proxy_row)
                    }

                    if (proxy_present == FALSE) {
                        print(paste0("No proxy SNP available for ", missing_IVs[i], " in outcome GWAS."))
                    }
                }
            }

            # Harmonising exposure and outcome datasets
            data_harmon <- harmonise_data(
                exposure_dat = exp_dat_clumped,
                outcome_dat = out_dat,
                action = 2
            )

            # Generate output file names
            harmon_data_file <- paste0("out/", "harmon_", file_name, sep = "")

            # Write the data to CSV files
            write.csv(data_harmon, harmon_data_file)
        },
        error = function(e) {
            # Handle the error
            error_file_name <- paste0("out/", "error_", file_name, sep = "")
            error_message <- paste0("Error: ", conditionMessage(e))
            error_data <- data.frame(File = file_name, Error = error_message)
            write.csv(error_data, error_file_name)
        }
    )
}
