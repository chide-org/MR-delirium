# 安装并加载必要的包
library(httr)
library(xml2)

# 要查询的 rs 号列表
exposure <- extract_instruments(outcomes = "ieu-b-110")
temp <- grepl("^rs*", dat$SNP) #正则表达式匹配SNP列中以rs开头的行
nrow(dat) # 总行数
sum(temp)
sum(!temp) # 不符合要求的总行数，大于0进行大括号里面的操作
exposure <- dat[temp]
rm(dat)
gc()
colnames(dat)
range = 1:50
computeREF(dat,rs_col = "SNP",ref_col = "A2",4000:5000,echo = F)
nrow(exposure)
sum(temp)
sum(!temp)

findREF <- function(rsid) {
  tryCatch({
    id <- substring(rsid, 3)
    id <- regmatches(id, gregexpr("[0-9]+", id))[[1]]
    # 构造查询的 URL
    base_url <-
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    response <- GET(
      url = paste0(
        base_url,
        "?db=snp&id=",
        id,
        "&api_key=fca0f05531dafacf5c1a1010c2c26ea93408",
        sep = ""
      ),
      timeout(5)
    )
    if (status_code(response) == 200) {
      # 解析返回的 XML 数据
      xml_content <- content(response, "parsed")
      parsed_xml <- read_xml(xml_content)
      
      # 提取 SNP 的参考位点信息
      parsed_xml <- xml_find_all(parsed_xml, ".//DOCSUM")
      if (length(parsed_xml) > 0) {
        ref <-
          regmatches(parsed_xml,
                     gregexpr("\\[[^\\[/]+/[^\\[]+\\]", parsed_xml))[[1]]
        return(ref)
      }
      else{
        print(paste(rsid, "not find in daSNP"))
        return("[?/?]")
      }
    } else {
      print(paste(rsid, ":Failed to fetch data, refetching..."))
      return("[*/*]")
    }
  },
  error = function(e) {
    print(paste0(rsid, e))
    return("[?/?]")
  })
}

computeREF <- function(data, rs_col="SNP", ref_col="A2", range=1:20) {
  # 确保 rs_col 和 ref_col 是字符向量
  rs_col <- as.character(rs_col)
  ref_col <- as.character(ref_col)
  
  #data <- data[grepl("^rs*", data[,..rs_col])]
  temp <- c()
  for (i in range) {
    SNP <- data[i, ..rs_col]
    if(!grepl("^rs*", SNP))next
    REF <- data[i, ..ref_col]
    text <- findREF(SNP)
    ref <- substring(text, 2, 2)
    if (ref == "*") {
      i <- i - 1
      next
    }
    if (ref == "?")
      next
    if (ref == REF) {
      bool <- T
    } else {
      bool <- F
      print(paste0(SNP, " : ", REF, "!=", ref, " ", bool))
    }
    temp <- c(temp, bool)
  }
  print(paste0("ref：", sum(temp)))
  print(paste0("non-ref：", sum(!temp)))
  print(paste0("referrence allele rate：", sum(temp) / (sum(temp) + sum(!temp))))
}
