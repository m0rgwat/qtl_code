arguments <- commandArgs(trailingOnly = TRUE)

library(DESeq2)
library(preprocessCore)

MedianNorm_EBSeq <- function (Data, alternative = FALSE) 
{
  if (ncol(Data) == 1) 
    stop("Only 1 sample!")
  if (!alternative) {
    geomeans <- exp(rowMeans(log(Data)))
    out <- apply(Data, 2, function(cnts) median((cnts/geomeans)[geomeans > 
                                                                  0]))
  }
  if (alternative) {
    DataMatO <- Data
    N <- ncol(DataMatO)
    DataList0 <- sapply(1:N, function(i) DataMatO[, i]/DataMatO, 
                        simplify = F)
    DataEachMed0 <- sapply(1:N, function(i) apply(DataList0[[i]], 
                                                  2, function(j) median(j[which(j > 0 & j < Inf)])))
    DataColgeo <- sapply(1:N, function(i) exp(mean(log(DataEachMed0[-i, 
                                                                    i]))))
    out <- DataColgeo
  }
  out
}

#also

quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

med_norm_Can <- function(data = data){
  
  row_med <- apply(data,1,function(x){median(x,na.rm = T)})
  med_div <- data/row_med
  size_factors <- apply(med_div,2,function(x){median(x[x>0],na.rm = T)})
  return(size_factors)
}

expression <- read.delim(file=arguments[1], sep = "\t", 
                         row.names = NULL, header = T,
                         stringsAsFactors = F)

sample_range <- c(7:ncol(expression))

expression[,sample_range] <- apply(expression[,sample_range],2,function(x){
  
  x <- as.numeric(x)
  
})

if (arguments[2]=="median_Can"){
  
  norm_factor <- med_norm_Can(expression[,sample_range])
  expression[,sample_range] <- as.matrix(sweep(expression[,sample_range],2,norm_factor,"/"))
  
} else if (arguments[2]=="median_Can_quantile"){
  
  norm_factor <- med_norm_Can(expression[,sample_range])
  expression[,sample_range] <- as.matrix(sweep(expression[,sample_range],2,norm_factor,"/"))
  expression_mat <- as.matrix(expression[,sample_range])
  expression[,sample_range] <- normalize.quantiles(expression_mat)
  
} else if (arguments[2]=="median_EBSeq"){
  
  norm_factor <- MedianNorm_EBSeq(expression[,sample_range])
  expression[,sample_range] <- as.matrix(sweep(expression[,sample_range],2,norm_factor,"/"))
  
} else if (arguments[2]=="median_EBSeq_quantile"){
  
  norm_factor <- MedianNorm_EBSeq(expression[,sample_range])
  expression[,sample_range] <- as.matrix(sweep(expression[,sample_range],2,norm_factor,"/"))
  expression_mat <- as.matrix(expression[,sample_range])
  expression[,sample_range] <- quantile_normalisation(expression_mat)
  
} else if (arguments[2]=="quantile"){
  
  expression_mat <- as.matrix(expression[,sample_range])
  expression[,sample_range] <- normalize.quantiles(expression_mat)
  
} else if (arguments[2]=="VST"){
  
  expression[,sample_range] <- varianceStabilizingTransformation(as.matrix(expression[,sample_range]),
                                                                 fitType = "mean")
  
} else if (arguments[2]=="none"){
  

} else if (arguments[2]=="median_Helene"){
  
  size_factors <- estimateSizeFactorsForMatrix(as.matrix(expression[,sample_range]), type = "poscounts")
  expression[,sample_range] <- as.matrix(sweep(expression[,sample_range],2,size_factors,"/"))
  
} else if (arguments[2]=="median_Helene_quantile"){
  
  size_factors <- estimateSizeFactorsForMatrix(as.matrix(expression[,sample_range]))
  expression_mat <- as.matrix(sweep(expression[,sample_range],2,size_factors,"/"))
  expression[,sample_range] <- quantile_normalisation(expression_mat)
  
}

colnames(expression)[1] <- "#chr"

write.table(expression, arguments[3], sep = "\t", quote = F, row.names = F)








