arguments <- commandArgs(trailingOnly = TRUE)

print("read nominal data")
nominal_data <- data.table::fread(arguments[1], sep = " ", 
                                  header = T, stringsAsFactors = F, data.table = F)

expression_data <- read.delim(arguments[2], sep = "\t", header = T)
n_genes <- nrow(expression_data)

print("Number of comparisons (n_genes*n_snips)")

print(paste0(n_genes, " * ", as.numeric(arguments[3]), " = ", n_genes*as.numeric(arguments[3])))

print("to numeric and remove NAs")
nominal_data$NPVAL <- as.numeric(nominal_data$NPVAL)
nominal_data <- nominal_data[!is.na(nominal_data$NPVAL),]

print("order pvalues")
nominal_data <- nominal_data[order(nominal_data$NPVAL),]

print("compute fdr")
nominal_data$fdr <- nominal_data$NPVAL/(c(1:nrow(nominal_data))/(n_genes*as.numeric(arguments[3])))

print("filter fdr <0.05")
nominal_data <- nominal_data[nominal_data$fdr<0.05,]

print(paste0("Number of eGenes (FDR <0.05): ", length(unique(as.character(nominal_data[,1])))))

print("writing the result")
write.table(nominal_data, arguments[4], quote=FALSE, sep = "\t", row.names = F)

