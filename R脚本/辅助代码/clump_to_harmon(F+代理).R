library(TwoSampleMR)
library(LDlinkR)
library(tidyverse)
library(dplyr)
library(stringr)
library(data.table)
# file_names_subset <- file_names
outcomePath <- "D:/A-MR/outcome/using_outcome.csv"
exposurePath <- "D:/A-MR/exposure/mydata/clump"
# 处理结局
outcomedat <- fread(outcomePath, header = TRUE)
viewn(outcomedat)
# 批量处理暴露
names <- list.files(exposurePath)

# Loop through each file
setwd(exposurePath)
for (file_name in names) {
    tryCatch(
        {
            # F
            EXP <- fread(file_name, header = TRUE)
            if (sum(is.na(EXP$eaf.exposure)) > 0) {
                print(paste0("eaf.exposure has NA in ", file_name))
                exp <- EXP
            } else {
                
                EXP$R2 <- (2 * EXP$beta.exposure * EXP$beta.exposure * EXP$eaf.exposure * (1 - EXP$eaf.exposure) / (2 * EXP$beta.exposure * EXP$beta.exposure * EXP$eaf.exposure * (1 - EXP$eaf.exposure) + 2 * EXP$se.exposure * EXP$se.exposure * EXP$samplesize.exposure * EXP$eaf.exposure * (1 - EXP$eaf.exposure))) # 计算R2
                EXP$F_value <- EXP$R2 * (EXP$samplesize.exposure - 2) / (1 - EXP$R2) # 计算F检验值
                exp <- EXP %>% filter(F_value > 10)
            }

            # proxy
            missing_IVs <- exp$SNP[!(exp$SNP %in% outcomedat$SNP)]


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
                        # proxy_row$beta.outcome <- as.numeric(outcomedat[outcomedat$SNP == proxy_SNP, "beta.outcome"])
                        # proxy_row$se.outcome <- as.numeric(outcomedat[outcomedat$SNP == proxy_SNP, "se.outcome"])
                        if (outcomedat[outcomedat$SNP == proxy_SNP, "effect_allele.outcome"] == proxy_SNP_allele_1) proxy_row$effect_allele.outcome <- original_SNP_allele_1
                        if (outcomedat[outcomedat$SNP == proxy_SNP, "effect_allele.outcome"] == proxy_SNP_allele_2) proxy_row$effect_allele.outcome <- original_SNP_allele_2
                        if (outcomedat[outcomedat$SNP == proxy_SNP, "other_allele.outcome"] == proxy_SNP_allele_1) proxy_row$other_allele.outcome <- original_SNP_allele_1
                        if (outcomedat[outcomedat$SNP == proxy_SNP, "other_allele.outcome"] == proxy_SNP_allele_2) proxy_row$other_allele.outcome <- original_SNP_allele_2
                        # proxy_row$pval.outcome <- as.numeric(outcomedat[outcomedat$SNP == proxy_SNP, "pval.outcome"])
                        # if ("eaf.outcome" %in% colnames(outcomedat)) proxy_row$eaf.outcome <- as.numeric(outcomedat[outcomedat$SNP == proxy_SNP, "eaf.outcome"])
                        outcomedat <- rbind(outcomedat, proxy_row)
                    }

                    if (proxy_present == FALSE) {
                        print(paste0("No proxy SNP available for ", missing_IVs[i], " in outcome GWAS."))
                    }
                }
            }

            # Harmonising exposure and outcome datasets
            data_harmon <- harmonise_data(
                exposure_dat = exp,
                outcome_dat = outcomedat,
                action = 2
            )

            # Generate output file names
            harmon_data_file <- paste0("harmon/", "h_", substring(file_name, 7), sep = "")

            # Write the data to CSV files
            fwrite(data_harmon, harmon_data_file, row.names = FALSE, nThread = 8)
        },
        error = function(e) {
            # Handle the error
            error_file_name <- paste0("harmon/", "error_", substring(file_name, 7), sep = "")
            error_message <- paste0("Error: ", conditionMessage(e))
            error_data <- data.frame(File = file_name, Error = error_message)
            write.csv(error_data, error_file_name)
        }
    )
}
