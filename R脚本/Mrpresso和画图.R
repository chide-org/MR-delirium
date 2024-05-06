
library(plyr)
library(TwoSampleMR)



#循环
file_names <- list.files(pattern = "^harmon.*\\.csv$")
file_sizes <- file.size(file_names)
file_names <- file_names[file_sizes > 3]

file_names_subset <- file_names

for (file_name in file_names_subset)  {
  tryCatch({
data_h <- read.csv(file_name)

presso <- run_mr_presso(data_h)

write.csv(presso[[1]]$`Main MR results`, paste0("mr_presso_Main", file_name))
write.csv(presso[[1]]$`MR-PRESSO results`$`Global Test`, paste0("mr_presso_Global", file_name))
write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, paste0("mr_presso_Outlier", file_name))

  }, error = function(e) {
    # Handle the error
    error_file_name <- paste0("error_", file_name)
    error_message <- paste0("Error: ", conditionMessage(e))
    error_data <- data.frame(File = file_name, Error = error_message)
    write.csv(error_data, error_file_name)
  })
}










#画图
dat <- read.csv("h_ebi-a-GCST006572.csv")


#绘制散点图
pdf(file="h_ebi-a-GCST006572.scatter_plot.pdf", width=7.5, height=7)
mrResult=mr(dat)
mr_scatter_plot(mrResult, dat)
dev.off()

#漏斗图
pdf(file="h_ebi-a-GCST006572.funnel_plot.pdf", width=7, height=6.5)
res_single=mr_singlesnp(dat)      #得到每个工具变量对结局的影响
mr_funnel_plot(singlesnp_results = res_single)
dev.off()

#留一法敏感性分析
pdf(file="h_ebi-a-GCST006572.leaveoneout.pdf", width=7, height=6.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
dev.off()




