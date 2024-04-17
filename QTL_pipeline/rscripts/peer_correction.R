arguments <- commandArgs(trailingOnly = TRUE)


expression <- read.delim(file=arguments[1], sep = "\t", 
                         row.names = NULL, header = T)

metadata <- read.csv(file=arguments[2], 
                     row.names = 1, header = T, stringsAsFactors = T)

sample_range <- c(7:ncol(expression))

expression_mat <- t(expression[,sample_range])


n_unobserved_factors <- round(nrow(expression_mat)*0.6,0)
if (n_unobserved_factors>100){n_unobserved_factors <- 100}


metadata[,] <- sapply(metadata,as.numeric)
metadata_mat <- as.matrix(metadata)

if (arguments[5]=="PEER"){
  
  library(peer)
  
  peer_model <- PEER()
  PEER_setNk(peer_model, n_unobserved_factors)
  PEER_setPhenoMean(peer_model, expression_mat)
  #PEER_setCovariates(peer_model, metadata_mat)
  
  PEER_update(peer_model)
  
  peer_factors <- PEER_getX(peer_model)
  
  rownames(peer_factors) <- rownames(expression_mat)
  colnames(peer_factors) <- paste0("peer_f_No_",1:ncol(peer_factors))
  
  pdf(file = paste0(arguments[4],".pdf"), width = 15, height = 15)
  PEER_plotModel(peer_model)
  dev.off()
  
} else if (arguments[5]=="PCA"){
  
  expression_mat <- scale(expression_mat,T,T)
  peer_factors <- prcomp(expression_mat, center = F)$x[,1:n_unobserved_factors]
  
  colnames(peer_factors) <- paste0("Expression_computed_",1:ncol(peer_factors))

}



peer_factors <- t(peer_factors)

peer_factors <- cbind(rep("chr1",n_unobserved_factors),
                      seq(10, by = 10, length.out = n_unobserved_factors),
                      seq(10, by = 10, length.out = n_unobserved_factors)+1,
                      rownames(peer_factors),
                      rownames(peer_factors),
                      rep("+", n_unobserved_factors),
                      peer_factors)

colnames(peer_factors)[1:6] <- c("#chr", "start", "end",
                                 "phenotype_id", "phenotype_group_id",
                                 "strand")


write.table(peer_factors,
            file = arguments[3],
            row.names = F, col.names = T, quote = F, sep = "\t")


