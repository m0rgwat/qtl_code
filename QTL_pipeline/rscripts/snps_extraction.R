arguments <- commandArgs(trailingOnly = TRUE)

conditional_results <- read.delim(file=arguments[1], sep = " ", 
                         row.names = NULL, header = F,
                         stringsAsFactors = F)

if (ncol(conditional_results)>=24&grepl(":",conditional_results[1,10], fixed = T)==T){
  
  snp_column <- 10
  
} else {
  
  snp_column <- 8
  
}

conditional_results <- conditional_results[!is.na(conditional_results[,snp_column]),]

snp_list <- do.call(rbind,sapply(conditional_results[,snp_column], function(x){
  
  strsplit(x,":",fixed = T)
  
}))

write.table(snp_list, arguments[2], sep = "\t", quote = F, row.names = F, col.names = F)








