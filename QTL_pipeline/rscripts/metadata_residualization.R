arguments <- commandArgs(trailingOnly = TRUE)


peer_ass_nominal <- read.delim(file=arguments[1], sep = " ", 
                               row.names = NULL, header = T)
peer_ass_permutational <- read.delim(file=arguments[2], sep = " ", 
                                     row.names = NULL, header = F)
expression <- read.delim(file=arguments[3], sep = "\t", 
                         row.names = NULL, header = T)
metadata <- read.csv(file=arguments[4], 
                     row.names = 1, header = T)
peer_factors <- read.delim(file=arguments[5], sep = "\t", 
                           row.names = NULL, header = T)
minimal_perm_p <- min(as.numeric(peer_ass_permutational[,3]))

peer_ass_nominal <- peer_ass_nominal[peer_ass_nominal$PCHR!="PCHR",]

peer_correlated <- peer_ass_nominal[as.numeric(peer_ass_nominal$NPVAL)<minimal_perm_p,]

if (nrow(peer_correlated)>0){
  peer_ass_nominal <- peer_ass_nominal[!(peer_ass_nominal[,1]%in%unique(peer_correlated[,1])),]}

uncorrelated_peer_factors <- peer_factors$phenotype_id[peer_factors$phenotype_id%in%unique(peer_ass_nominal[,1])]

write.table(peer_correlated, 
            paste0(arguments[5],"_ass_w_genotype"), sep = "\t", quote = F, row.names = T)

write.table(peer_factors$phenotype_id[!(peer_factors$phenotype_id%in%uncorrelated_peer_factors)], 
            paste0(arguments[5],"_ass_names"), sep = "\t", quote = F, row.names = T)


peer_factors_unc <- t(peer_factors[,7:ncol(peer_factors)])
colnames(peer_factors_unc) <- peer_factors$phenotype_id
rownames(peer_factors_unc) <- colnames(peer_factors)[7:ncol(peer_factors)]
peer_factors <- peer_factors_unc


peer_factors <- peer_factors[,uncorrelated_peer_factors]

metadata <- metadata[,-1]

if (arguments[7]=="auto") {
  
  if (isTRUE(all(rownames(metadata)==rownames(peer_factors)))){
    
    metadata <- cbind(metadata,peer_factors)
    
  } else {stop("Metadata and peer factors have different rownames")}
  
} else 
  
  if (as.numeric(arguments[7])>0&as.numeric(arguments[7])<ncol(peer_factors)) {
  
  peer_factors <- peer_factors[,1:arguments[7]]
  
  if (isTRUE(all(rownames(metadata)==rownames(peer_factors)))){
    
    metadata <- cbind(metadata,peer_factors)
    
  } else {stop("Metadata and peer factors have different rownames")}
  
} else 
  
  if (as.numeric(arguments[7])==0) {} else {
    
    if (isTRUE(all(rownames(metadata)==rownames(peer_factors)))){
      
      metadata <- cbind(metadata,peer_factors)
      
    } else {stop("Metadata and peer factors have different rownames")}
    
}




sample_range <- c(7:ncol(expression))

expression[,sample_range] <- t(scale(apply(t(expression[,sample_range]), 2, function(x){
  
  residuals(lm(x~., data = metadata))
  
})))

colnames(expression)[1] <- "#chr"

write.table(expression, arguments[6], sep = "\t", quote = F, row.names = F)



