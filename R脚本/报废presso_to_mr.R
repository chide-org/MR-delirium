library(tidyverse)
library(formattable)
library(htmltools)
library(webshot)
library(qpdf) 
export_formattable <- function(f, file, width = "100%", height = NULL,
                               background = "white", delay = 0.2) {
    w <- as.htmlwidget(f, width = width, height = height)
    path <- html_print(w, background = background, viewer = NULL)
    url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
    webshot(url,
        file = file,
        selector = ".formattable_widget",
        delay = delay
    )
}

# formatter
mr1 <- formatter("span",
    style = ~ ifelse(pval > 0.05, "color:red", "color:green")
) ## 大于0.05显示黑色，否则绿色
mr2 <- formatter("span",
    style = ~ ifelse(b > 0, "color:green", "color:red")
)
mr3 <- formatter("span",
    style = ~ ifelse(method == "Inverse variance weighted", "color:red", "color:black")
)
h <- formatter("span",
    style = ~ ifelse(Q_pval < 0.05, "color:red", "color:black")
)
p <- formatter("span",
    style = ~ ifelse(pval < 0.05, "color:red", "color:black")
)

# prepare data
harmonDirs <- "exposure/harmon/"
files <- list.files(harmonDirs, full.names = F)
files

# cycle deal data
for (file in files) {
    dat_harmon <- fread(paste0(harmonDirs, file))
    file <- substring(file, 8)
    name <- substring(file, 0, nchar(file) - 4)

    # mr
    mr <- mr(dat_harmon)
    hetero <- mr_heterogeneity(dat_harmon)
    pleio <- mr_pleiotropy_test(dat_harmon)

    mr_table <- formattable(mr, list(pval = mr1, b = mr2, method = mr3))
    hetero_table <- formattable(hetero, list(Q_pval = h))
    pleio_table <- formattable(pleio, list(pval = p))

    export_formattable(mr_table, paste0("mroutput/", name, "_mr", ".pdf", sep = ""))
    export_formattable(hetero_table, paste0("mroutput/", name, "_hetero", ".pdf", sep = ""))
    export_formattable(pleio_table, paste0("mroutput/", name, "_pleio", ".pdf", sep = ""))
    # plot
    pdf(file = paste0("mroutput/", name, ".pdf", sep = ""), height = 5, width = 8)
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
    print(p1)
    print(p2)
    print(p3)
    dev.off()
    pdf_combine(c(paste0("mroutput/", name, "_mr", ".pdf", sep = ""),paste0("mroutput/", name, "_hetero", ".pdf", sep = ""),
                  paste0("mroutput/", name, "_pleio", ".pdf", sep = ""),paste0("mroutput/", name, ".pdf", sep = "")),
                output = paste0("mroutput/", name, ".pdf", sep = ""))
    file.remove(paste0("mroutput/", name, "_mr", ".pdf", sep = ""),paste0("mroutput/", name, "_hetero", ".pdf", sep = ""))
    file.remove(paste0("mroutput/", name, "_pleio", ".pdf", sep = ""))   
}
