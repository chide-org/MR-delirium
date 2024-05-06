library(tidyverse)


# prepare data
setwd("D:/A-MR/exposure/harmon")
harmonDirs <- "/"
files <- list.files(harmonDirs, full.names = F)
files



# cycle deal data
for (file in files) {
    tryCatch(
        {
            dat_harmon <- fread(file)
            file <- substring(file, 8)
            name <- substring(file, 0, nchar(file) - 4)


            # mr
            mr <- mr(dat_harmon)
            write.csv(mr, file = paste0("../out/", "mr_", file, sep = ""))

            hetero <- mr_heterogeneity(dat_harmon)
            write.csv(hetero, file = paste0("../out/", "hetero_", file, sep = ""))

            pleio <- mr_pleiotropy_test(dat_harmon)
            write.csv(pleio, file = paste0("../out/", "pleio_", file, sep = ""))
            # plot

            p1 <- mr_scatter_plot(
                mr_results = mr(
                    dat_harmon,
                    method_list = c(
                        "mr_ivw",
                        "mr_weighted_median",
                        "mr_egger_regression"
                    )
                ), dat_harmon
            )
            p2 <- mr_funnel_plot(singlesnp_results = mr_singlesnp(dat_harmon))
            p3 <- mr_leaveoneout_plot(
                leaveoneout_results =
                    mr_leaveoneout(dat_harmon)
            )
            pdf(file = paste0("../out/", "egger_", name, ".pdf", sep = ""), height = 8, width = 12)
            print(p1)
            dev.off()

            pdf(file = paste0("../out/", "funnel_", name, ".pdf", sep = ""), height = 8, width = 12)
            print(p2)
            dev.off()

            pdf(file = paste0("../out/", "leaveoneout_", name, ".pdf", sep = ""), height = 8, width = 12)
            print(p3)
            dev.off()
        },
        error = function(e) {
            print(paste0("File: ", file))
            print(paste0("Error: ", e))
        }
    )
}
