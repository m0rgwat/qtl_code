arguments <- commandArgs(trailingOnly = TRUE)

library(vcfR)

expression <- read.delim(file=arguments[2], sep = "\t", 
                         row.names = NULL, header = T)

sample_range <- c(7:ncol(expression))
gene_names <- 4

cis_on_step <- data.table::fread(file = arguments[1],
                                 sep = " ", header = FALSE,
                                 select = c(1, 8), 
                                 col.names = c("tag_ID", "var_ID"))

if (grepl("chr",unlist(cis_on_step[1,2]))){}else{
  
  cis_on_step <- data.table::fread(file = arguments[1],
                                   sep = " ", header = FALSE,
                                   select = c(1, 10), 
                                   col.names = c("tag_ID", "var_ID"))
  
}

vfc_012 <- read.vcfR(file = arguments[3])

################################################################################
##############################_vcfR2dosages_####################################

vcfR2dosages <- function (x, return.alleles = F) 
{
  x <- extract.gt(x, return.alleles = return.alleles, element = "DS")
  x <- as.data.frame(t(x))
  icol <- 1:ncol(x)
  for (i in icol) x[, i] <- factor(x[, i])
  class(x) <- c("loci", "data.frame")
  attr(x, "locicol") <- icol
  x
}


################################################################################

vfc_012 <- vcfR2dosages(vfc_012, return.alleles = F)
vfc_012 <- apply(vfc_012,c(1,2),function(x){as.numeric(as.character(x))})

expression_correction <- as.data.frame(t(expression[,sample_range]))
vfc_012 <- vfc_012[match(row.names(expression_correction),row.names(vfc_012)),]
colnames(expression_correction) <- expression[,gene_names]

for (i in 1:nrow(cis_on_step)){
  
  expr_on_step <- expression_correction[,colnames(expression_correction)==cis_on_step$tag_ID[i]]
  
  if (is.vector(expr_on_step)==T){
    
    expression_correction[,colnames(expression_correction)==cis_on_step$tag_ID[i]] <- resid(
      lm(expr_on_step~vfc_012[,colnames(vfc_012)==cis_on_step$var_ID[i]]))
    
  } else {
    
    expression_correction[,colnames(expression_correction)==cis_on_step$tag_ID[i]] <- apply(
      expr_on_step,2, function(x){
        
        resid(lm(x~vfc_012[,colnames(vfc_012)==cis_on_step$var_ID[i]]))
      })
    
  }
  

  
}

expression[,sample_range] <- t(scale(expression_correction, T, T))
colnames(expression)[1] <- "#chr"

write.table(expression,
            file = arguments[4],
            row.names = F, col.names = T, quote = F, sep = "\t")

