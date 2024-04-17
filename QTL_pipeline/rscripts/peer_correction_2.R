arguments <- commandArgs(trailingOnly = TRUE)

library(peer)

expression <- read.delim(file=arguments[1], sep = "\t", 
                         row.names = NULL, header = T)

metadata <- read.csv(file=arguments[2], 
                     row.names = 1, header = T, stringsAsFactors = F)

sample_range <- c(7:ncol(expression))

expression_mat <- t(expression[,sample_range])


n_unobserved_factors <- round(nrow(expression_mat)*0.3,0)
if (n_unobserved_factors>100){n_unobserved_factors <- 100}


metadata[,] <- apply(metadata,2,as.numeric)
metadata_mat <- as.matrix(metadata)

peer_model <- PEER()
PEER_setNk(peer_model, n_unobserved_factors)
PEER_setPhenoMean(peer_model, expression_mat)
PEER_setCovariates(peer_model, metadata_mat)
PEER_setNmax_iterations(peer_model, as.numeric(arguments[4]))

PEER_update(peer_model)

peer_residuals <- PEER_getResiduals(peer_model)

expression[,sample_range] <- t(peer_residuals)

pdf(file = paste0("peer_model.pdf"), width = 15, height = 15)
PEER_plotModel(peer_model)
dev.off()

colnames(expression)[1] <- "#chr"
write.table(expression,
            file = arguments[3],
            row.names = F, col.names = T, quote = F, sep = "\t")


