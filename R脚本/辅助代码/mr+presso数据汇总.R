setwd("D:/A-MR/exposure")
presso <- read.csv("无混淆_presso_output.csv",header = T)

mr <- read.csv("无混淆_mr_output.csv",header = T)
all <- left_join(mr,presso,by = "ID.exposure")

write.csv(all,"无混淆_mr_output.csv",row.names = F)
