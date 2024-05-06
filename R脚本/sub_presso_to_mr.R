library(tidyverse)
library(gridExtra)

# prepare data
setwd("D:/A-MR/exposure/remr")
harmonDirs <- "./"
files <- list.files(harmonDirs, full.names = F)
files


# cycle deal data
for (file in files) {
    tryCatch(
        {
            dat_harmon <- fread(file)
            file <- substring(file, 9)
            name <- substring(file, 0, nchar(file) - 4)


            # mr
            mr <- mr(dat_harmon)
            hetero <- mr_heterogeneity(dat_harmon)
            pleio <- mr_pleiotropy_test(dat_harmon)
            # plot
            pdf(file = paste0("../out/", name, ".pdf", sep = ""), height = 8, width = 12)

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
            empty_row4 <- data.frame(matrix("", nrow = 20, ncol = ncol(mr)))
            colnames(empty_row4) <- colnames(mr)
            mr <- mr %>%
                rbind(empty_row4) %>%
                select(5:ncol(mr))
            p4 <- grid.table(mr)
            print(p4)

            if (nrow(hetero) != 0) {
                hetero <- hetero %>%
                    mutate(type = "heterogeneity", .before = "method")
                empty_row5 <- data.frame(matrix("", nrow = 8, ncol = ncol(hetero)))
                colnames(empty_row5) <- colnames(hetero)
                hetero <- hetero %>%
                    rbind(empty_row5) %>%
                    select(5:ncol(hetero))
                p5 <- grid.table(hetero)
                print(p5)
            }

            if (nrow(pleio) != 0) {
                pleio <- pleio %>%
                    select(4:ncol(pleio)) %>%
                    mutate(type = "pleiotropy_test") %>%
                    select(type, everything())
                p6 <- grid.table(pleio)
                print(p6)
            }
            print(p1)
            print(p2)
            print(p3)

            dev.off()
        },
        error = function(e) {
            print(paste0("File: ", file))
            print(paste0("Error: ", e))
        }
    )
}
