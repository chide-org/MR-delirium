# 教程：https://mp.weixin.qq.com/s/F9BA_YI-RG-IomNvBKKLTw
# download Rtools
# https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html

# download SNPlocs.Hsapiens.dbSNP155.GRCh37和SNPlocs.Hsapiens.dbSNP155.GRCh38
# http://bioconductor.org/packages/release/data/annotation/html/SNPlocs.Hsapiens.dbSNP155.GRCh37.html
# http://bioconductor.org/packages/release/data/annotation/html/SNPlocs.Hsapiens.dbSNP155.GRCh38.html

SNP_LOC_DATA <- SNPlocs.Hsapiens.dbSNP155.GRCh38::SNPlocs.Hsapiens.dbSNP155.GRCh38

# 一键运行
{
  partdat <- dat[temp]
  wrongdat <- dat[!temp, ] %>% select(-SNP)
  wrongdat$CHR <- as.character(wrongdat$CHR)
  wrongdat$BP <- as.character(wrongdat$BP)
  wrongdat <- wrongdat[complete.cases(wrongdat[, c("CHR", "BP"), with = FALSE])]
  wrongdat$CHR <- ifelse(wrongdat$CHR %in% c("23", "X", "x"), "X", wrongdat$CHR)
  wrongdat <- split(wrongdat, wrongdat$CHR) # 按照染色体分组
  gc()
  file_name <- c()
  setwd(exposurePath)
  getwd()
  for (wdat in wrongdat) {
    file_path <- paste0("chr_rs/", wdat[1, "CHR"], ".csv")
    data.table::fwrite(wdat, file = file_path, row.names = FALSE)
    file_name <- c(file_name, file_path)
  }

  # 主要步骤，耗时
  sapply(file_name, function(file) {
    wdat <- data.table::fread(file, data.table = F)
    rsids <- wdat %>%
      GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE, seqnames.field = "CHR", start.field = "BP", end.field = "BP")
    rsids <- data.table::setDT(data.frame(BSgenome::snpsByOverlaps(SNP_LOC_DATA, ranges = rsids)))
    rsids <- rsids %>%
      dplyr::filter(!duplicated(paste0(rsids$seqnames, rsids$pos))) %>%
      dplyr::select(seqnames, pos, RefSNP_id) %>%
      dplyr::rename(CHR = seqnames, BP = pos, SNP = RefSNP_id)
    rsids$CHR <- as.character(rsids$CHR)
    rsids$BP <- as.character(rsids$BP)
    wdat$CHR <- as.character(wdat$CHR)
    wdat$BP <- as.character(wdat$BP)
    dplyr::left_join(wdat, rsids, by = c("CHR", "BP")) %>% data.table::fwrite(file, row.names = F)
    gc()
    return(NULL)
  })

  # 合并
  wdat <- lapply(file_name, function(file) {
    wdat <- data.table::fread(file)
  }) %>% data.table::rbindlist()

  file.remove(file_name)

  # 去除空列
  subtemp <- grepl("^rs*", wdat$SNP)
  rightdat <- wdat[subtemp]

  nrow(partdat)
  nrow(rightdat)

  dat <- rbind(partdat, rightdat) %>% unique()
  nrow(dat)
}
