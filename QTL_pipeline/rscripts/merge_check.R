arguments <- commandArgs(trailingOnly = TRUE)

expression <- read.delim(file=arguments[1], sep = "\t", 
                         row.names = NULL, header = T)

colnames(expression)[1] <- "#chr"

metadata <- read.csv(file=arguments[2],  
                     row.names = NULL, header = T, stringsAsFactors = F)

#metadata <- metadata[metadata[,2]==arguments[3],]
rownames(metadata) <- metadata[,1]

sample_list <- colnames(expression)[7:ncol(expression)]

metadata <- metadata[rownames(metadata)%in%sample_list,]
metadata <- metadata[match(sample_list,row.names(metadata)),]

if (all(rownames(metadata)==sample_list)){"Everything ok"} else {
  
  print("Can not find sample(s): ")
  print(rownames(metadata)[!(rownames(metadata)%in%sample_list)])
  print("in the metadata file")
  
  stop(print("Row names in metadata and expression data differs!"),call. = T)
  
}

if ((colnames(expression)[6]=="strand"&
    all(colnames(expression)[1:5]==c("#chr","start",
                                 "end","phenotype_id",
                                 "phenotype_group_id")))){"Everything ok"} else {
  
  print("There is a problem in the header of .BED file: ")
  print(arguments[1])
  
  stop(print("Please check the header!"),call. = T)
  
}

if (grepl("chr", unlist(expression[1,1]))){"Everything ok"} else {
  
  print("Your chromosomes names are incorrect, you have: ")
  print(unlist(expression[1,1]))
  print("But it should starts with chr!")
  print("The problem is found for: ")
  print(arguments[1])
  stop(print("Please change the chromosomes' names!"),call. = T)
                                       
 }

write.csv(metadata, paste0("metadata_merged_", arguments[3]), quote = F, row.names = F)
write.table(sample_list,"sample_list", sep = "\t", quote = F, row.names = F, 
            col.names = F)



