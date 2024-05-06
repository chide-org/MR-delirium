rm(list = ls())
gc()
library(dplyr)
library(data.table)
library(stringr)
library(data.table)
library(magrittr)
# download Rtools
# https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html

# download SNPlocs.Hsapiens.dbSNP155.GRCh37和SNPlocs.Hsapiens.dbSNP155.GRCh38
# http://bioconductor.org/packages/release/data/annotation/html/SNPlocs.Hsapiens.dbSNP155.GRCh37.html
# http://bioconductor.org/packages/release/data/annotation/html/SNPlocs.Hsapiens.dbSNP155.GRCh38.html

SNP_LOC_DATA <- SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38
# 获取数据
data <- fread("outcome/outcome_empty.csv") %>%
    dplyr::rename(
        CHR = `#chrom`, # 对应染色体
        BP = pos
    ) # 对应pos
gc()
data$CHR <- as.character(data$CHR)
data$BP <- as.character(data$BP)
gc()
data <- data[complete.cases(data[, c("CHR", "BP"), with = FALSE])]
gc()
data$CHR <- ifelse(data$CHR %in% c("23", "X", "x"), "X", data$CHR)
gc()
data <- split(data, data$CHR)#按照染色体分组
gc()
file_name <- c()
sapply(data, function(dat) {
    data.table::fwrite(dat, file = paste0("temp/",nrow(dat), ".csv"), row.names = F)
    file_name <- c(file_name, paste0("temp/", nrow(dat), ".csv"))
})#将分组数据分别保存到本地

sapply(file_name, function(file) {
    data <- data.table::fread(file, data.table = F)
    rsids <- data %>%
        GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE, seqnames.field = "CHR", start.field = "BP", end.field = "BP")
    rsids <- data.table::setDT(data.frame(BSgenome::snpsByOverlaps(SNP_LOC_DATA, ranges = rsids)))
    gc()
    rsids <- rsids %>%
        dplyr::filter(!duplicated(paste0(rsids$seqnames, rsids$pos))) %>%
        dplyr::select(seqnames, pos, RefSNP_id) %>%
        dplyr::rename(CHR = seqnames, BP = pos, SNP = RefSNP_id)
    rsids$CHR <- as.character(rsids$CHR)
    rsids$BP <- as.character(rsids$BP)
    data$CHR <- as.character(data$CHR)
    data$BP <- as.character(data$BP)
    dplyr::left_join(data, rsids, by = c("CHR", "BP")) %T>% data.table::fwrite(file, row.names = F)
    gc()
    return(NULL)
})

data <- lapply(file_name, function(file) {
    data <- data.table::fread(file)
}) %>% data.table::rbindlist()

file.remove(file_name)


viewn(data)
View(data)
