############_LOAD_PACKAGES_#################
# setwd("/home/vpetrov/r_work/theta_cor/theta_function")

library(Rfast)
library(ggplot2)
library(scales)
library(data.table)
library(vcfR)
library(ggpubr)
library(biomaRt)
library(rtracklayer)


############################################

#########_START_DEFINGING_FUNCTIONS_########

############_FAST_LM_FUNCTION_##############

# CHECKED AND FIXED

fast_lm <- function (y = y, x = x) {
  
  n <- dim(x)[1]
  denom <- n - 1
  my <- sum(y)/n
  m <- Rfast::colmeans(x)
  r <- (Rfast::eachcol.apply(x, y) - n * my * m)/denom
  sx <- Rfast::colVars(x, suma = n * m)
  be <- r/sx
  bex <- t(t(x) * be)
  a <- my - be*m
  pred_y <- bex + rep(a, each = nrow(bex))
  
  Se <- sqrt((Rfast::colsums((t(t(x)*be+a)-y)^2)/(n-2))/Rfast::colsums((t(t(x)-m))^2))
  
  pvalue <- pt(abs(be/Se), n-2, lower.tail = FALSE)*2
  result <- cbind(be, pvalue)
  
  if (is.null(colnames(x))) {
    rownames(result) <- paste("X", 1:ncol(x), sep = "")
  }
  else rownames(result) <- colnames(x)
  
  return(result)
  
}

############################################

############_FAST_COR_FUNCTION_#############

# CHECKED AND FIXED

fastcor <- function (y, x, type = "pearson") {
  
  n <- length(y)
  if (type == "pearson") {
    r <- as.vector(cor(y, x))
  }
  else if (type == "spearman") {
    r <- as.vector(cor(Rfast::Rank(y), Rfast::colRanks(x)))
  }
  tval <- r/sqrt((1-r^2)/(n-2))
  pvalue <- 2 * pt(abs(tval), n - 2, lower.tail = FALSE)
  res <- cbind(r, pvalue)
  colnames(res) <- c("correlation", "p-value")
  if (is.null(colnames(x))) {
    rownames(res) <- paste("X", 1:ncol(x), sep = "")
  }
  else rownames(res) <- colnames(x)
  res
  
}

############################################

###############_DATA_LOAD_##################

# CHECKED AND FIXED

load_data <- function (file = file, locus = locus) {
  
  loaded_data <- readRDS(file)
  
  if (is.data.frame(loaded_data)){} else {
    
    loaded_data <- loaded_data[[locus]]
    
  }
  return(loaded_data)
  
}

############################################

########_FAST_CORRELATION_OR_LM_############

# CHECKED AND FIXED

correl_on_step <- function(snp_mat = snp_mat, 
                           expr_vec = expr_vec, 
                           cortype = cortype){
  
  ## snp_mat - matrix with SNPs dosages, samples per row, SNPs per column
  ## expr_vec - vector with expression values for a gene inside locus
  ## cortype - type of model used. Supported types are "lm", "pearson" and "spearman"
  
  
  if (cortype=="lm") {
    
    res <- fast_lm(y = expr_vec, x = snp_mat)
    
    colnames(res) <- c("cor_or_beta","p_value")
    return(res)
    
  } else {
    
    res <- fastcor(y = expr_vec, x = snp_mat,
                          type = cortype)
    colnames(res) <- c("cor_or_beta","p_value")
    return(res)
    
  }
  
}

############################################


############_THETA_CORRELATION_#############

# CHECKED AND FIXED

theta_cor <- function(X_p = X_p, X_beta = NULL, 
                      Y_p = Y_p, Y_beta = NULL, 
                      useSign_beta = useSign) {
  
  if (isTRUE(useSign_beta)){
    
    X_logp <- -log10(X_p)
    Y_logp <- -log10(Y_p)

  } else {
    
    X_logp <- -log10(X_p)
    Y_logp <- -log10(Y_p)
    
    X_beta <- 1
    Y_beta <- 1
    
  }
  
  
  w_i <- mapply(max, X_logp/max(X_logp), Y_logp/max(Y_logp))
  
  X_w <- sum(X_logp*w_i)/sum(w_i)
  Y_w <- sum(Y_logp*w_i)/sum(w_i)
  
  X_sigma <- sqrt(sum(w_i*(X_logp-X_w)^2)/sum(w_i))
  Y_sigma <- sqrt(sum(w_i*(Y_logp-Y_w)^2)/sum(w_i))
  
  R_w <- sum(w_i*((X_logp-X_w)/X_sigma)*((Y_logp-Y_w)/Y_sigma))/sum(w_i)
  R_ws <- sum(w_i*((X_logp-X_w)/X_sigma)*((Y_logp-Y_w)/Y_sigma)*sign(X_beta*Y_beta))/sum(w_i)
  
  theta <- R_ws/((1+exp(1)^(-30*(R_w-0.3))))
  
  return(theta)
  
}


#  more fast (to check)

theta_cor <- function(X_p = X_p, X_beta = NULL, 
                      Y_p = Y_p, Y_beta = NULL, 
                      useSign_beta = useSign) {
  
  X_logp <- -log10(X_p)
  Y_logp <- -log10(Y_p)
  
  if (isTRUE(!useSign_beta)){
    
    X_beta <- 1
    Y_beta <- 1
    
  }
  
  
  w_i <- pmax(X_logp/max(X_logp), Y_logp/max(Y_logp))
  
  X_w <- weighted.mean(X_logp,w_i)
  Y_w <- weighted.mean(Y_logp,w_i)
  
  ## weighted.var modi - to replace
  X_sigma <- sqrt(sum(w_i*(X_logp-X_w)^2)/sum(w_i))
  Y_sigma <- sqrt(sum(w_i*(Y_logp-Y_w)^2)/sum(w_i))
  
  ## weighted.cor weights
  R_w <- sum(w_i*((X_logp-X_w)/X_sigma)*((Y_logp-Y_w)/Y_sigma))/sum(w_i)
  ## the same, but multiply by sign
  R_ws <- sum(w_i*((X_logp-X_w)/X_sigma)*((Y_logp-Y_w)/Y_sigma)*sign(X_beta*Y_beta))/sum(w_i)
  
  theta <- R_ws/(1+exp(-30*(R_w-0.3)))
  
  return(theta)
  
}

############################################

###############_CHECK_DATA_#################

#########_CHECK_EXPRESSION_GENOTYPE_########

# CHECKED AND FIXED

check_eg <- function(expression_data = expression_data,
                     locus_genotype = locus_genotype){
  
  col_intersect <- intersect(colnames(expression_data[,-1]),
                             colnames(locus_genotype[,3:ncol(locus_genotype)]))
  
  expression_data <- expression_data[,c("PID",col_intersect)]
  locus_genotype <- locus_genotype[,c("SNPname","ID",col_intersect)]
  
  results <- list(expression_data,locus_genotype)
  names(results) <- c("expression_data","locus_genotype")
  
  
  return(results)
  
}

############################################

#########_CHECK_GENOTYPE_GWAS_##############

# CHECKED AND FIXED

check_gg <- function(DA = DA,
                     locus_genotype = locus_genotype){
  
  
  intersection <- intersect(DA[,"SNPname"],locus_genotype[,"SNPname"])
  
  DA <- DA[DA$SNPname%in%intersection,]
  DA <- DA[match(intersection,DA$SNPname),]
  
  locus_genotype <- locus_genotype[locus_genotype$SNPname%in%intersection,]
  locus_genotype <- locus_genotype[match(intersection,locus_genotype$SNPname),]
  
  results <- list(DA,locus_genotype)
  names(results) <- c("DA","locus_genotype")
  
  
  return(results)
  
}

############################################

############_CHECK_TWO_GWASes_##############

# CHECKED AND FIXED

check_gwgw <- function(DA1 = DA1,
                       DA2 = DA2){
  
  
  intersection <- intersect(DA1[,"SNPname"],DA2[,"SNPname"])
  
  DA1 <- DA1[DA1$SNPname%in%intersection,]
  DA1 <- DA1[match(intersection,DA1$SNPname),]
  
  DA2 <- DA2[DA2$SNPname%in%intersection,]
  DA2 <- DA2[match(intersection,DA2$SNPname),]
  
  results <- list(DA1,DA2)
  names(results) <- c("DA1","DA2")
  
  
  return(results)
  
}

############################################

#########_CIS_CIS_INTERSECTIONS_############

# CHECKED AND FIXED

cis_cis_intersections <- function(cis_results = cis_results,
                                  sym_border = sym_border){
  
  cis_results_test <- cis_results[order(cis_results$p_value),]
  cis_results_test$tg_identifier <- paste0(cis_results_test$PID,":",cis_results_test$tissue)
  results <- as.data.frame(matrix(nrow = 0,
                                  ncol = 9))
  
  for (i in 1:nrow(cis_results)){
    
    chr_test <- cis_results_test$snp_chr[i]
    adapt_border <- c(cis_results_test$snp_pos[i]-sym_border, 
                      cis_results_test$snp_pos[i]+sym_border)
    
    if (adapt_border[1]<0){adapt_border[1] <- 1}
    
    cis_subset <- cis_results_test[cis_results_test$snp_chr==chr_test&
                                     cis_results_test$snp_pos>=adapt_border[1]&
                                     cis_results_test$snp_pos<=adapt_border[2]&
                                     cis_results_test$tg_identifier!=cis_results_test$tg_identifier[i],]
    
    if (nrow(cis_subset)>0){
      
      left_border <- sapply(cis_subset$snp_pos, function(x){min(c(adapt_border,max(c(x-sym_border,1)),(x+sym_border)))})
      right_border <- sapply(cis_subset$snp_pos, function(x){max(c(adapt_border,max(c(x-sym_border,1)),(x+sym_border)))})
      
      ## test data.table::rbindlist()
      
      results <- rbind(results,list(rep(cis_results_test$PID[i],nrow(cis_subset)),
                                    cis_subset$PID,
                                    rep(cis_results_test$snp_ID[i],nrow(cis_subset)),
                                    cis_subset$snp_ID,
                                    cis_subset$snp_chr,
                                    left_border,right_border,
                                    rep(cis_results_test$tissue[i],nrow(cis_subset)),
                                    cis_subset$tissue))
    }
    
  }
  
  colnames(results) <- c("target_PID", "recip_PID",
                         "target_snp_ID", "recip_snp_ID", "chr", 
                         "left_border","right_border","target_tissue","recip_tissue")
  
  for_duplicates <- results[,c(1,2,8,9)]
  for_duplicates <- t(apply(for_duplicates,1,function(x){
    
    c(paste0(x[1],"_",x[3]),paste0(x[2],"_",x[4]))
    
  }))
  
  duplicates <- apply(for_duplicates,1,function(x){
    
    paste(c(unlist(x)[order(unlist(x))]),collapse = "")
    
  })
  
  results <- results[!duplicated(duplicates),]
  
  return(results)
  
}

############################################

#########_CIS_CIS_INTERSECTIONS_############

DAP_mapping <- function(loci = loci,
                        cis_results = cis_results,
                        sym_border = sym_border){
  
  results <- as.data.frame(matrix(nrow = 0,ncol = 9))
  
  for (i in 1:nrow(loci)){
    
    genes_inside <- cis_results[cis_results$snp_chr==loci$chr[i]&
                                  cis_results$snp_pos>=loci$left_border[i]&
                                  cis_results$snp_pos<=loci$right_border[i],]
    
    
    
    if (nrow(genes_inside)>0){
      
      left_border <- sapply(genes_inside$snp_pos, function(x){min(c(genes_inside$left_border,genes_inside$right_border,
                                                                    max(c(x-sym_border,1)),(x+sym_border)))})
      left_border_min <- min(left_border)
      
      right_border <- sapply(genes_inside$snp_pos, function(x){max(c(genes_inside$left_border,genes_inside$right_border,
                                                                     max(c(x-sym_border,1)),(x+sym_border)))})
      right_border_max <- max(right_border)
      
      ## test data.table::rbindlist()
      
      results <- rbind(results, list(rep(loci$rsID[i],nrow(genes_inside)),
                                     genes_inside$PID,
                                     genes_inside$snp_ID,
                                     genes_inside$p_value,
                                     rep(loci$chr[i],nrow(genes_inside)),
                                     rep(left_border_min,nrow(genes_inside)),
                                     rep(right_border_max,nrow(genes_inside)),
                                     left_border,
                                     right_border))
      
    }
    
  }
  
  colnames(results) <- c("locus_rsID", "PID",
                         "cis_snp_ID", "p_value", "chr", "left_border","right_border", 
                         "left_gene_specific","right_gene_specific")
  
  return(results)
  
}

############################################

###########_LOAD_NOMINAL_DATA_##############

# CHECKED AND FIXED

load_nominal <- function(nominal_data = nominal_data,
                         mapping = mapping,
                         analysis = analysis){
  
  if (analysis == "eap_eap"){
    
    mapping_transf <- mapping[,c(1,5:8)]
    colnames(mapping_transf)[c(1,5)] <- colnames(mapping)[c(2,9)]
    mapping_transf <- rbind(mapping_transf,
                            mapping[,c(2,5:7,9)])
    
    duplicates <- apply(mapping_transf,1,function(x){
      
      paste(gsub(" ", "", c(unlist(x))),collapse = "")
      
    })
    
    mapping_transf <- mapping_transf[!duplicated(duplicates),]
    mapping_transf <- mapping_transf[mapping_transf$recip_tissue%in%unique(nominal_data$tissue),]
    
    EAP <- apply(mapping_transf,1,function(x){
      
      as.data.frame(nominal_data[nominal_data$snp_chr==x[2]&
                                   nominal_data$snp_pos>=as.numeric(x[3])&
                                   nominal_data$snp_pos<=as.numeric(x[4])&
                                   nominal_data$PID==x[1]][,c("ID","cor_or_beta","p_value")])
      
    }, simplify = list)
    
    names(EAP) <- apply(mapping_transf, 1, function(x){
      
      gsub(" ", "", paste(x, collapse = ":"))
      
    })
    
    EAP <- lapply(EAP, function(x){
      
      colnames(x)[1] <- "ID"
      x$SNPname <- sapply(x$ID,ID_to_SNPname)
      x
      
    })
    
  } else if (analysis == "dap") {
    
    mapping_transf <- mapping[,c(5:7)]
    duplicates <- apply(mapping_transf,1,function(x){
      
      paste(gsub(" ", "", c(unlist(x))),collapse = "")
      
    })
    
    mapping_transf <- mapping_transf[!duplicated(duplicates),]
    
    EAP <- apply(mapping_transf,1,function(x){
      
      as.data.frame(nominal_data[nominal_data$chr==x[1]&
                                   nominal_data$pos>=as.numeric(x[2])&
                                   nominal_data$pos<=as.numeric(x[3])][,c("ID","SNPname",
                                                                              "cor_or_beta",
                                                                              "p_value")])
      
    }, simplify = list)
    
    names(EAP) <- apply(mapping_transf, 1, function(x){
      
      gsub(" ", "", paste(x, collapse = ":"))
      
    })

  }
  

  return(EAP)
  
}

############################################

#########_LOAD_EXPRESSION_DATA_#############

# CHECKED AND FIXED

load_expression <- function(expression_data = expression_data,
                            mapping = mapping, 
                            analysis = analysis){
  
  ## analysis - "eap_eap"; "dap"
  
  if (analysis == "eap_eap"){
    
    tg_pairs_target <- apply(mapping,1,function(x){
      
      paste0(x[1],":",x[8])
      
    })
    
    tg_pairs_recip <- apply(mapping,1,function(x){
      
      paste0(x[2],":",x[9])
      
    })
    
    tg_pairs <- c(tg_pairs_target,tg_pairs_recip)
    
    gene_expression <- list()
    
    for (tiss in unique(unlist(mapping[,c(8:9)]))){
      
      if ((!is.data.frame(expression_data))&is.list(expression_data)){
        
        expression_data_step <- expression_data[[tiss]]
        
      } else {expression_data_step <- expression_data}
      
      
      tg_pairs_sel <- sapply(tg_pairs,function(x){
        
        grepl(tiss,strsplit(x,":")[[1]][2])
        
      })
      
      tg_pairs_step <- tg_pairs[tg_pairs_sel]
      
      gene_expression_step <- sapply(unique(tg_pairs_step), function(x){
        
        y <- strsplit(x,":")[[1]][1]
        
        on_step <- expression_data_step[expression_data_step$phenotype_group_id==y, 
                                        c(5,7:ncol(expression_data_step)), drop = F]
        colnames(on_step)[1] <- "PID"
        on_step
        
      }, simplify = F)
      
      names(gene_expression_step) <- unique(tg_pairs_step)
      gene_expression <- append(gene_expression,gene_expression_step)
      
    }
    
  }
  
  if (analysis == "dap"){
    
    gene_expression <- sapply(unique(mapping$locus_rsID), function(x){
      
      on_step <- expression_data[expression_data$phenotype_group_id%in%
                                   mapping[mapping$locus_rsID%in%x,"PID"], 
                                      c(5,7:ncol(expression_data)), drop = F]
      colnames(on_step)[1] <- "PID"
      on_step
      
    }, simplify = F)
    
  }
  
  return(gene_expression)
  
}

############################################

#########_LOAD_GENOTYPES_DATA_#############

# CHECKED AND FIXED

load_genotypes <- function(vcf_file = vcf_file, mapping = mapping){
  
  ## analysis - "eap_eap"; "dap"
  
  genotypes <- read.vcfR(vcf_file)
  dosages <- vcfR2dosages(genotypes)
  
  gen_coordinates <- as.data.frame(t(sapply(colnames(dosages), function(x){
    
    y <-c(x, strsplit(x,":")[[1]][1], strsplit(x,":")[[1]][2])
    y
    
  })))
  
  colnames(gen_coordinates) <- c("ID","chr","pos")
  gen_coordinates$pos <- as.numeric(gen_coordinates$pos)
  
  mapping_locus <- mapping[,c(5:7)]
  
  duplicates <- apply(mapping_locus,1,function(x){
    
    paste(gsub(" ", "", c(unlist(x))),collapse = "")
    
  })
  
  mapping_locus <- mapping_locus[!duplicated(duplicates),]
  
  loci <- apply(mapping_locus,1,function(x){
    
    locus <- t(dosages[,gen_coordinates$chr==x[1]&
                         gen_coordinates$pos>=as.numeric(x[2])&
                         gen_coordinates$pos<=as.numeric(x[3]), drop = F])
    
    ID <- rownames(locus)
    locus <- as.data.frame(apply(locus, 2, as.numeric))
    locus <- cbind(ID,locus)
    
  }, simplify = F)
  
  names(loci) <- apply(mapping_locus, 1, function(x){
    
    gsub(" ", "", paste(x, collapse = ":"))
    
  })
  
  
  loci <- lapply(loci, function(x){
    
    x$SNPname <- sapply(x$ID,ID_to_SNPname)
    x <- x[,c(1,ncol(x),2:(ncol(x)-1))]
    x
    
  })

  
  return(loci)
  
}

############################################


###########_MAPPING_FOR_VCF_################

# CHECKED AND FIXED

make_VCF_map <- function(mapping = mapping){
  
  VCF_map <- mapping[order(sapply(mapping$chr, 
                                  function(x){as.numeric(sub("chr","",sub("chrX","23",x)))}),
                           mapping$left_border),]
  
  VCF_map <- VCF_map[,c("chr","left_border","right_border")]
  
  include_reg <- rep(T, nrow(VCF_map))
  
  for (i in 1:(nrow(VCF_map)-1)){
    
    if (VCF_map[i,"right_border"]>VCF_map[i+1,"left_border"]&
        VCF_map[i,"chr"]==VCF_map[i+1,"chr"]){
      
      VCF_map[i+1,"left_border"] <- VCF_map[i,"left_border"]
      
      include_reg[i] <- F
    }
    
  }
  
  VCF_map <- VCF_map[include_reg,]
  return(VCF_map)
  
}

############################################

#########_CHECK_REF_ALT_ALLELE_#############

# CHECKED AND FIXED

ref_alt_check <- function(locus_genotype = locus_genotype,
                          DA = DA){
  
  for (i in 1:nrow(locus_genotype)){
    
    if (tolower(strsplit(locus_genotype[i,"ID"],":")[[1]][[3]])==tolower(strsplit(DA[i,"ID"],":")[[1]][[3]])){} else {
      
      DA[i,"cor_or_beta"] <- DA[i,"cor_or_beta"]*(-1)
      
      }
  }
  return(DA)
}

############################################

#########_FILTER_WEAK_ASSOCIATIONS_#########

# CHECKED AND FIXED

f_weak_ass <- function(locus_genotype = locus_genotype,
                      expression_data = expression_data,
                      p_threshold = p_threshold,
                      cortype = cortype){
  
  ## locus_genotype - data.frame of genotypes dosages 
  ## with patients in columns and SNPs in rows.
  ## First two columns must be SNPname and ID
  ## Order of patients in columns should be identical 
  ## with the order in gene_expression
  
  ## expression_data - data.frame (even with 1 row)
  ## with the firs column named PID and containing gene names;
  ## genes expression are in rows, patients are in columns.
  
  ## p_threshold - borderline value of -log10(p) of the best association 
  ## to decide to keep the gene or to remove
  
  ## cortype - type of correlation used to compute expression
  ## association pattern (EAP) 
  ## choose from "lm", "spearman", "pearson"
 
  locus_genotype_val <- t(locus_genotype[,-c(1,2)])
  p_threshold <- -log10(p_threshold)
  min_genes_p <- sapply(expression_data$PID, function(x){
    
    min(correl_on_step(snp_mat = locus_genotype_val,
                       expr_vec = as.numeric(as.character(unlist(expression_data[expression_data$PID==x,-1]))),
                       cortype = cortype)[,2], na.rm = T)
    
  })
  
  if (all((-log10(min_genes_p))<p_threshold)){
    
    return(print(paste0("No genes found with -log10(p) >= ", p_threshold)))
    
  } else {
    
    expression_data <- expression_data[expression_data$PID%in%names(min_genes_p[(-log10(min_genes_p))>=p_threshold]),, drop = F]
    
  }
  
  return(expression_data)
  
}

############################################

#########_NOMINAL_P_FOR_THETA_##############

# CHECKED AND FIXED

step_theta <- function(locus_genotype = locus_genotype,
                       expression_data = expression_data,
                       cortype = cortype, DA = DA,
                       useSign = useSign,
                       return_p_vec = return_p_vec){
  
  locus_genotype_val <- t(locus_genotype[,-c(1,2)])
  
  colnames(locus_genotype_val) <- locus_genotype$SNPname
  
  correl_res <- list()
  
  correl_res <- sapply(expression_data$PID, function(x){
    
    on_step <- list(correl_on_step(snp_mat = locus_genotype_val,
                   expr_vec = as.numeric(as.character(unlist(expression_data[expression_data$PID==x,-1]))),
                   cortype = cortype))
    
    correl_res <- append(correl_res,on_step)

  })
  
  
  nominal_theta_R <- sapply(expression_data$PID, function(x){
    
    correl_on_step <- correl_res[[x]]
    
    theta_cor(X_p = DA$p_value, X_beta = DA$cor_or_beta, 
              Y_p = correl_on_step[,"p_value"], Y_beta = correl_on_step[,"cor_or_beta"],
              useSign = useSign)
    
  })
  
  if (return_p_vec == T){
    
    return(list(nominal_theta_R,correl_res))
    
  } else {
    
    return(nominal_theta_R)
    
    }
}

############################################

########_PERMUTATIONAL_ENGINE_##############

# CHECKED AND FIXED

permutheta <- function(locus_genotype = locus_genotype,
                       expression_data = expression_data,
                       cortype = cortype, DA = DA,
                       B = B, useSign = useSign){
  
  nominal_theta_R <- step_theta(locus_genotype = locus_genotype,
                                expression_data = expression_data,
                                cortype = cortype, DA = DA,
                                useSign = useSign,
                                return_p_vec = F)
  
  theta_permutational_p <- rep(1,length(nominal_theta_R))
  names(theta_permutational_p) <- names(nominal_theta_R)

  for (i in c(1:B)){
    
    expression_data_perm <- expression_data[,c("PID",
                                          sample(colnames(expression_data[-1]),
                                                          replace = F)), drop = F]
    
    theta_permutational_step <- step_theta(locus_genotype = locus_genotype,
                                  expression_data = expression_data_perm,
                                  cortype = cortype, DA = DA,
                                  useSign = useSign,
                                  return_p_vec = F)

    
    theta_permutational_p <- theta_permutational_p + 
      (abs(theta_permutational_step) > abs(nominal_theta_R))
    
  }
  
  theta_permutational_p <- theta_permutational_p/(B+1)
  
  result <- cbind(nominal_theta_R,theta_permutational_p)
  colnames(result)[1:2] <- c("theta_R","p")
  
  return(result)
  
}

############################################

########_FILTER_UNSIGNIFICANT_SNPs_#########

filt_SNPs <- function(locus_genotype = locus_genotype,
                      DA = DA,
                      nominal_res = nominal_res,
                      SNP_threshold = SNP_threshold){
  
  if (is.null(dim(nominal_res[,grepl("p_value",colnames(nominal_res))]))){
    
    genes_p_include <- nominal_res[,grepl("p_value",colnames(nominal_res))]<=SNP_threshold
    
  } else {
    
    genes_p_include <- apply(nominal_res[,grepl("p_value",colnames(nominal_res))],1,function(x){
      
      any(unlist(x)<=SNP_threshold)
      
    })
    
  }
  
  SNPs_to_keep <- genes_p_include|(DA$p_value<=SNP_threshold)
  
  if (sum(SNPs_to_keep)<10){
    
    stop(print(paste0("Analysis stopped! 
                        To few SNPs with p-value less than ", SNP_threshold)))
    
  }
  
  locus_genotype <- locus_genotype[SNPs_to_keep,]
  DA <- DA[SNPs_to_keep,]
  
  return(list(locus_genotype,DA))
}

############################################

######_ADAPTIVE_PERMUTATIONAL_ENGINE_#######

# CHECKED AND FIXED

permutheta_adaptive <- function(locus_genotype = locus_genotype,
                                expression_data = expression_data,
                                cortype = cortype, DA = DA,
                                B = B, useSign = useSign, 
                                adaptive = adaptive, alpha = alpha,
                                SNP_threshold = SNP_threshold){
  
  ## locus_genotype - data.frame of genotypes dosages 
  ## with patients in columns and SNPs in rows.
  ## First two columns must be SNPname and ID
  ## Order of patients in columns should be identical 
  ## with the order in gene_expression
  
  ## expression_data - data.frame (even with 1 row)
  ## with the firs column named PID and containing gene names;
  ## genes expression are in rows, patients are in columns.
  
  ## cortype - type of correlation used to compute expression
  ## association pattern (EAP) 
  ## choose from "lm", "spearman", "pearson"
  
  ## DA - GWAS results EAP or DAP
  
  ## B - number of permutations
  
  ## useSign - usage of correlation sign in theta computation
  
  ## adaptive usage of adaptive permutation procedure
  
  ## alpha - significance border after association is treated as
  ## non-significant
  
  ## SNP_threshold <- threshold for SNPs to be included in theta calculation
  
  nominal_theta <- as.data.frame(step_theta(locus_genotype = locus_genotype,
                                expression_data = expression_data,
                                cortype = cortype, DA = DA,
                                useSign = useSign,
                                return_p_vec = T)[[2]])
  
  SNPs_filtered <- filt_SNPs(locus_genotype = locus_genotype,
                             DA = DA,
                             nominal_res = nominal_theta,
                             SNP_threshold = SNP_threshold)
  
  locus_genotype <- SNPs_filtered[1][[1]]
  DA <- SNPs_filtered[2][[1]]
  
  nominal_theta_R <- step_theta(locus_genotype = locus_genotype,
                                     expression_data = expression_data,
                                     cortype = cortype, DA = DA,
                                     useSign = useSign,
                                     return_p_vec = F)
  
  theta_permutational_p <- matrix(rbind(rep(1,length(nominal_theta_R)), 
                                        rep(0,length(nominal_theta_R))),
                                  ncol = length(nominal_theta_R), nrow = 2)
  
  colnames(theta_permutational_p) <- names(nominal_theta_R)
  row.names(theta_permutational_p) <- c("counts","b")
  
  early_stop <- round((B+1)*alpha, digits = 0)+1
  
  b <- 0
  
  while (b!=B&!(nrow(expression_data)==0&adaptive)){
    
    expression_data_perm <- expression_data[,c("PID",
                                               sample(colnames(expression_data[-1]),
                                                      replace = F)), drop = F]
    
    theta_permutational_step <- step_theta(locus_genotype = locus_genotype,
                                           expression_data = expression_data_perm,
                                           cortype = cortype, DA = DA,
                                           useSign = useSign,
                                           return_p_vec = F)
    
    theta_permutational_p[1,expression_data$PID] <- theta_permutational_p[1,expression_data$PID] + 
      (abs(theta_permutational_step) > abs(nominal_theta_R[expression_data$PID]))
    
    theta_permutational_p[2,expression_data$PID] <- theta_permutational_p[2,expression_data$PID]+1
    
    
    if (any(theta_permutational_p[1,]>=early_stop)){
      
      to_exclude <- colnames(theta_permutational_p)[which(theta_permutational_p[1,]>=early_stop)]
      expression_data <- expression_data[!(expression_data$PID%in%to_exclude), ,drop = F]
      
    }
    
    b <- b+1
    
  }
  
  theta_permutational_p[1,] <- theta_permutational_p[1,]/(theta_permutational_p[2,]+1)
  
  result <- cbind(nominal_theta_R,t(theta_permutational_p))
  colnames(result)[1:3] <- c("theta_R","p","b")
  
  return(result)
  
}

############################################

###########_ID_TO_SNPname_##################

# CHECKED AND FIXED

ID_to_SNPname <- function(x){
  
  name_ends <- strsplit(x,":")[[1]][3:4]
  name_ends <- name_ends[order(name_ends)]
  
  return(paste(c(strsplit(x,":")[[1]][1:2], 
                 name_ends), collapse = ":"))
  
}

############################################

###########_VCF_TO_DOSAGES_#################

# CHECKED AND FIXED

vcfR2dosages <- function (x, return.alleles = F) 
{
  x <- extract.gt(x, return.alleles = return.alleles, element = "DS")
  x <- as.data.frame(t(x))
  icol <- 1:ncol(x)
  class(x) <- c("loci", "data.frame")
  attr(x, "locicol") <- icol
  x[,] <- apply(x,2,function(y){as.numeric(y)})
  x
}

############################################

########_THETA_CORRELATION_EAP_EAP_#########

# CHECKED AND FIXED

theta_eap_eap <- function(data_target = data_target,
                          data_recip = data_recip,
                          locus_genotype = locus_genotype,
                          cortype = cortype,
                          B = B, useSign = useSign, 
                          adaptive = adaptive, 
                          alpha = alpha,
                          SNP_threshold){
  
  ## locus_genotype - data.frame of genotypes dosages 
  ## with patients in columns and SNPs in rows.
  ## First two columns must be SNPname and ID
  ## Order of patients in columns should be identical 
  ## with the order in gene_expression
  
  ## data_target - list with data for "target" gene 
  ## (simply the first gene in mapping). It contains 
  ## two leafs named as follow:
  
  ## expression_data - data.frame (even with 1 row)
  ## with the firs column named PID and containing gene names;
  ## genes expression are in rows, patients are in columns.
  #
  ## DA - eQTL result (EAP) for this gene
  
  ## data_recip - list with data for "recipient" gene 
  ## (simply the second gene in mapping). It contains 
  ## two leafs named as follow:
  
  ## expression_data - data.frame (even with 1 row)
  ## with the firs column named PID and containing gene names;
  ## genes expression are in rows, patients are in columns.
  #
  ## DA - eQTL result (EAP) for this gene
  
  ## cortype - type of correlation used to compute expression
  ## association pattern (EAP) 
  ## choose from "lm", "spearman", "pearson"
  
  ## B - number of permutations
  
  ## useSign - usage of correlation sign in theta computation
  
  ## adaptive usage of adaptive permutation procedure
  
  ## alpha - significance border after association is treated as
  ## non-significant
  
  ## SNP_threshold <- threshold for SNPs to be included in theta calculation
  
  ###################################################
  
  ## prepare the data
  
  gwgw_checked <- check_gwgw(DA1 = data_target[["DA"]], 
                             DA2 = data_recip[["DA"]])
  
  expression_data_r <- data_recip[["expression_data"]]
  EAP_t <- gwgw_checked[[1]]
  locus_genotype_r <- locus_genotype
  
  eg_checked_r <- check_eg(expression_data = expression_data_r,
                           locus_genotype = locus_genotype_r)
  
  locus_genotype_r <- eg_checked_r[[2]]
  expression_data_r <- eg_checked_r[[1]]
  
  
  gg_checked_t <- check_gg(DA = EAP_t, 
                           locus_genotype = locus_genotype_r)
  
  EAP_t <- gg_checked_t[[1]]
  locus_genotype_r <- gg_checked_t[[2]]
  

  
  expression_data_t <- data_target[["expression_data"]]
  EAP_r <- gwgw_checked[[2]]
  locus_genotype_t <- locus_genotype
  
  eg_checked_t <- check_eg(expression_data = expression_data_t,
                           locus_genotype = locus_genotype_t)
  
  locus_genotype_t <- eg_checked_t[[2]]
  expression_data_t <- eg_checked_t[[1]]
  
  
  gg_checked_r <- check_gg(DA = EAP_r, 
                           locus_genotype = locus_genotype_t)
  
  EAP_r <- gg_checked_r[[1]]
  locus_genotype_t <- gg_checked_r[[2]]
  
  ###################################################
  
  ## run the analysis
  
  # permuting recip gene, so receive result for precomputed target
  res_target <- permutheta_adaptive(locus_genotype = locus_genotype_r,
                                    expression_data = expression_data_r,
                                    cortype = cortype, DA = EAP_t,
                                    B = B, useSign = useSign, 
                                    adaptive = adaptive, alpha = alpha,
                                    SNP_threshold = SNP_threshold)
  
  # permuting target gene, so receive result for precomputed recip
  res_recip <- permutheta_adaptive(locus_genotype = locus_genotype_t,
                                   expression_data = expression_data_t,
                                   cortype = cortype, DA = EAP_r,
                                   B = B, useSign = useSign, 
                                   adaptive = adaptive, alpha = alpha,
                                   SNP_threshold = SNP_threshold)
  
  ## end
  
  ###################################################
  
  result <- as.data.frame(matrix(ncol = 10, nrow = 1))
  colnames(result) <- c("Gene_target","Gene_recip",
                        "R_target","R_recip",
                        "b_target","b_recip",
                        "p_target","p_recip",
                        "R", "p_value")
  
  result[1,] <- list(rownames(res_recip),
                     rownames(res_target),
                     res_target[,1],res_recip[,1],
                     res_target[,3],res_recip[,3],
                     res_target[,2],res_recip[,2],
                     mean(c(res_target[,1],res_recip[,1])),
                     mean(c(res_target[,2],res_recip[,2])))
  
  return(result)
  
}


############################################

########_THETA_CORRELATION_EAP_EAP_#########


theta_dap <- function(expression_data = expression_data,
                      locus_genotype = locus_genotype,
                      DAP = DAP,
                      cortype = cortype,
                      B = B, useSign = useSign, 
                      adaptive = adaptive, 
                      alpha = alpha, p_threshold = p_threshold,
                      SNP_threshold = SNP_threshold){
  
  ## locus_genotype - data.frame of genotypes dosages 
  ## with patients in columns and SNPs in rows.
  ## First two columns must be SNPname and ID
  ## Order of patients in columns should be identical 
  ## with the order in gene_expression
  
  ## expression_data - data.frame (even with 1 row)
  ## with the firs column named PID and containing gene names;
  ## genes expression are in rows, patients are in columns.
  
  ## DAP - GWAS result (DAP) for the trait
  
  ## cortype - type of correlation used to compute expression
  ## association pattern (EAP) 
  ## choose from "lm", "spearman", "pearson"
  
  ## B - number of permutations
  
  ## useSign - usage of correlation sign in theta computation
  
  ## adaptive usage of adaptive permutation procedure
  
  ## alpha - significance border after association is treated as
  ## non-significant
  
  ## SNP_threshold <- threshold for SNPs to be included in theta calculation
  
  ###################################################
  
  ## prepare the data
  
  eg_checked <- check_eg(expression_data = expression_data,
                         locus_genotype = locus_genotype)
  
  expression_data <- eg_checked[[1]]
  locus_genotype <- eg_checked[[2]]
  
  gg_checked <- check_gg(DA = DAP, 
                         locus_genotype = locus_genotype)
  
  DAP <- gg_checked[[1]]
  locus_genotype <- gg_checked[[2]]
  
  DAP <- ref_alt_check(locus_genotype = locus_genotype,
                       DA = DAP)
  
  expression_data <- f_weak_ass(locus_genotype = locus_genotype,
                                expression_data = expression_data,
                                p_threshold = p_threshold,
                                cortype = cortype)
  
  if (typeof(expression_data)=="character"){
    
    return(print("Analysis stopped"))
    
  }
  
  ###################################################
  
  result <- permutheta_adaptive(locus_genotype = locus_genotype,
                                expression_data = expression_data,
                                cortype = cortype, DA = DAP,
                                B = B, useSign = useSign, 
                                adaptive = adaptive, alpha = alpha,
                                SNP_threshold = SNP_threshold)
  
  ###################################################
  
  return(result)
  
}



###################################################

###############_MAKE_TRAINING_LOCUS_###############

training_locus <- function(GWAS_dat = GWAS_dat,
                           chr = chr,
                           left_border = left_border,
                           right_border = right_border,
                           padding = padding){
  
  padding_size <- (right_border-left_border)*padding
  
  GWAS_subset <- as.data.frame(GWAS_dat[GWAS_dat$chr==chr&
                                          GWAS_dat$pos>=(left_border-padding_size)&
                                          GWAS_dat$pos<=(right_border+padding_size),])
  
  GWAS_subset$landscape <- 0
  GWAS_subset[(GWAS_subset$pos>=left_border)&(GWAS_subset$pos<=right_border),"landscape"] <- 1
  
  GWAS_subset <- GWAS_subset[,c("ID","SNPname","landscape","p_value","pos","chr")]
  return(GWAS_subset)
  
}

###################################################


#################_PRINTING_FUNCTIONS_##############

############_PRINTING_OF_SINGLE_GWAS_PEAK_#########

# CHECKED AND FIXED

print_gwas_peak <- function(GWAS_dat = GWAS_dat,
                            chr = chr,
                            snp_pos = snp_pos,
                            sym_border = sym_border,
                            adapt_border = NULL,
                            snp_ID = NULL,
                            draw_lines = NULL){
  
  ## GWAS_dat - file with GWAS summary statistics,
  ## should contain at leas these columns:
  ## p_value, pos, chr
  
  ## chr - chromosome of top SNP ("chr5") or NULL
  ## snp_pos - position of top SNP or NULL
  
  ## snp_ID - instead of chr and snp_pos you can 
  ## provide snp_ID in this format: chr5:56068166:C:T
  ## or "chr":"snp_pos":"not really interesting"
  ## You can't have both chr, snp_pos and snp_ID
  ## NULL or equal to something, one should be NULL
  
  ## sym_border - border of region to visualize;
  ## it is symmetric border around the position
  ## of top SNP. Should provide number of bases
  ## +/- from the top SNP
  ## (region would be sym_border*2) or NULL
  
  ## adapt_border - border of region to visualize;
  ## it is adaptive border, you can provide different 
  ## values for the left and the right borders. Could be
  ## c(left,right) in exact coordinates on chromosome 
  ## in number of bases or NULL
  
  ## You can't have both adapt_border and sym_border
  ## NULL or equal to something, one should be NULL
  
  ## draw_lines - optionally you can draw vertical
  ## lines. You should provide x-coordinate of 
  ## lines in bp as c(pos_1,pos_2,..) or NULL
  
  if (is.null(sym_border)&is.null(adapt_border)){
    
    stop("You should provide or symmetric or adaptive borders")
    
  } else if ((!is.null(sym_border))&(!is.null(adapt_border))) {
    
    stop("You can't provide both symmetric and adaptive borders")
    
  }
  
  if ((is.null(chr)|is.null(snp_pos))&is.null(snp_ID)){
    
    stop("You should provide or chr, snp_pos or snp_ID")
    
  } else if (((!is.null(chr))|(!is.null(snp_pos)))&(!is.null(snp_ID))){
    
    stop("You can't provide both chr, snp_pos and snp_ID")
    
  }
  
  if (!is.null(snp_ID)){
    
    chr <- strsplit(snp_ID,":")[[1]][1]
    snp_pos <- as.numeric(strsplit(snp_ID,":")[[1]][2])
    
  }
  
  if (is.null(adapt_border)){
    
    adapt_border <- c(snp_pos-sym_border, snp_pos+sym_border)
    
  }
  
  GWAS_subset <- as.data.frame(GWAS_dat[GWAS_dat$chr==chr&
                                          GWAS_dat$pos>=adapt_border[1]&
                                          GWAS_dat$pos<=adapt_border[2],])
  
  GWAS_subset$pos <- GWAS_subset$pos/1000
  
  gwas_peak <- ggplot(GWAS_subset, aes(x = pos, y = -log10(p_value))) +
    geom_point() + theme_bw() + 
    labs(title = paste0("GWAS peak for ", chr ," and position ", snp_pos))+
    xlab("Position, kb") + ylab("-log10(p-value)") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 15),
                       labels = scales::number_format(accuracy = 1, decimal.mark = '.')) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 11)) +
    theme(axis.text.x = element_text(angle = 15, size=12),
          axis.text.y = element_text(size=18)) + geom_hline(yintercept = (-log10(5e-8)))
  
  if (!is.null(draw_lines)){
    
    draw_lines <- draw_lines/1000
    gwas_peak <- gwas_peak + geom_vline(xintercept = draw_lines, 
                                        color = "blue", linetype = "dashed")
    
  }
  
  plot(gwas_peak)
  
  
}

##################################################


###############_PRINTING_MOUNTAINS_################

####################_THETA_PLOT_###################


plot_theta <- function(EAP_1 = EAP_1,
                       EAP_2 = EAP_2,
                       expression_data_1 = expression_data_1,
                       expression_data_2 = expression_data_2,
                       locus_genotype = locus_genotype,
                       cortype = cortype,
                       useSign = useSign,
                       gene_positions = gene_positions,
                       plot_type = plot_type,
                       draw_lines = draw_lines,
                       gene_1 = gene_1, 
                       gene_2 = gene_2,
                       tissue_1 = tissue_1, 
                       tissue_2 = tissue_2,
                       two_axes = two_axes,
                       weight_cor = weight_cor){
  
  
  ###################################################
  ### PLEASE IF YOU PROVIDE EAP 1 USE EXPRESSION DATA 2 SLOT
  ### AND VICE VERSA !!!!!!!!!!!!!!
  ###################################################
  
  ## EAP 1 and 2 - expression (disease) associated patterns with
  ## "gene" columns along with standard (SNPname, ID cor_or_beta, p_value). 
  ## Gene column contains the gene name
  
  ## expression_data 1 and 2 - data.frame with 1 row
  ## with the first column named PID and containing gene name;
  ## gene expression is in row, patients are in columns.
  
  # if one or both "expression_data" provided:
  #
  ## locus_genotype - data.frame of genotypes dosages 
  ## with patients in columns and SNPs in rows.
  ## First two columns must be SNPname and ID
  ## Order of patients in columns should be identical 
  ## with the order in gene_expression
  #
  ## cortype - type of correlation used to compute expression
  ## association pattern (EAP) 
  ## choose from "lm", "spearman", "pearson"
  #
  ## useSign - usage of correlation sign in theta computation
  
  ## draw_lines - optionally you can draw vertical
  ## lines. You should provide x-coordinate of 
  ## lines in bp as c(pos_1,pos_2,..) or NULL
  
  
  ## gene_positions - gene positions in current build
  
  ## plot_type - type of plot: "range" or "correlation"
  
  ## tissue_names - tissue names for EAP1 and EAP2
  ## gene - gene names for EAP1 and EAP2
  
  ## weight_cor - make plot with weighted values (T) or not (F)
  
  ###################################################
  
  if ((!is.null(expression_data_1)|!is.null(expression_data_2))&any(is.null(locus_genotype),
                                                                    is.null(cortype),
                                                                    is.null(useSign))){
    
    stop("If expression data is provided 
    you must provide locus_genotype, cortype, and useSign")
    
  }
  
  ###################################################
  
  if (((!is.null(EAP_1)&!is.null(EAP_2))&any(!is.null(expression_data_1),
                                             !is.null(expression_data_2)))|
      ((!is.null(EAP_1)&!is.null(EAP_2))&all(is.null(expression_data_1),
                                             is.null(expression_data_2)))){
    
    print("EAP_1 and EAP_2 will be used to make theta plot")
    
    ##############_EQUALIZE_EAP1_AND_EAP2_#############
    
    gwgw_checked <- check_gwgw(DA1 = EAP_1, DA2 = EAP_2)
    
    ###################################################
    
  } else if (is.null(EAP_1)&is.null(EAP_2)&
             (!is.null(expression_data_1)&!is.null(expression_data_2))){
    
    print("expression data will be used to make theta plot")
    
    ####################_MAKE_EAP1_####################
    
    locus_genotype_1 <- locus_genotype
    
    eg_checked_1 <- check_eg(expression_data = expression_data_1,
                             locus_genotype = locus_genotype_1)
    
    expression_data_1 <- eg_checked_1[[1]]
    locus_genotype_1 <- eg_checked_1[[2]]
    
    locus_genotype_val_1 <- t(locus_genotype_1[,-c(1,2)])
    colnames(locus_genotype_val_1) <- locus_genotype_1$SNPname
    
    EAP_1 <- as.data.frame(correl_on_step(snp_mat = locus_genotype_val_1,
                                          expr_vec = as.numeric(as.character(unlist(expression_data_1[,-1]))),
                                          cortype = cortype))
    
    EAP_1$SNPname <- locus_genotype_1$SNPname
    EAP_1$ID <- locus_genotype_1$ID
    
    ####################_MAKE_EAP2_####################
    
    locus_genotype_2 <- locus_genotype
    
    eg_checked_2 <- check_eg(expression_data = expression_data_2,
                             locus_genotype = locus_genotype_2)
    
    expression_data_2 <- eg_checked_2[[1]]
    locus_genotype_2 <- eg_checked_2[[2]]
    
    locus_genotype_val_2 <- t(locus_genotype_2[,-c(1,2)])
    colnames(locus_genotype_val_2) <- locus_genotype_2$SNPname
    
    EAP_2 <- as.data.frame(correl_on_step(snp_mat = locus_genotype_val_2,
                                          expr_vec = as.numeric(as.character(unlist(expression_data_2[,-1]))),
                                          cortype = cortype))
    
    EAP_2$SNPname <- locus_genotype_2$SNPname
    EAP_2$ID <- locus_genotype_2$ID
    
    ##############_EQUALIZE_EAP1_AND_EAP2_#############
    
    gwgw_checked <- check_gwgw(DA1 = EAP_1, DA2 = EAP_2)
    
    ###################################################
    
  } else if ((!is.null(EAP_1)|!is.null(EAP_2))&
             (!is.null(expression_data_1)|!is.null(expression_data_2))){
    
    print("expression data and EAP combination will be used to make theta plot")
    
    if (!is.null(EAP_1)&!is.null(expression_data_1)|
        !is.null(EAP_2)&!is.null(expression_data_2)){
      
      stop("Spirits of Machine are angry! PLEASE IF YOU PROVIDE EAP 1 USE EXPRESSION DATA 2 SLOT AND VICE VERSA !!!!!!!!!!!!!!")
      
    }
    
    if (!is.null(expression_data_1)){
      
      ####################_MAKE_EAP1_####################
      
      locus_genotype_1 <- locus_genotype
      
      eg_checked_1 <- check_eg(expression_data = expression_data_1,
                               locus_genotype = locus_genotype_1)
      
      expression_data_1 <- eg_checked_1[[1]]
      locus_genotype_1 <- eg_checked_1[[2]]
      
      locus_genotype_val_1 <- t(locus_genotype_1[,-c(1,2)])
      colnames(locus_genotype_val_1) <- locus_genotype_1$SNPname
      
      EAP_1 <- as.data.frame(correl_on_step(snp_mat = locus_genotype_val_1,
                                            expr_vec = as.numeric(as.character(unlist(expression_data_1[,-1]))),
                                            cortype = cortype))
      
      EAP_1$SNPname <- locus_genotype_1$SNPname
      EAP_1$ID <- locus_genotype_1$ID
      
      ###################################################
      
    } else {
      
      ####################_MAKE_EAP2_####################
      
      locus_genotype_2 <- locus_genotype
      
      eg_checked_2 <- check_eg(expression_data = expression_data_2,
                               locus_genotype = locus_genotype_2)
      
      expression_data_2 <- eg_checked_2[[1]]
      locus_genotype_2 <- eg_checked_2[[2]]
      
      locus_genotype_val_2 <- t(locus_genotype_2[,-c(1,2)])
      colnames(locus_genotype_val_2) <- locus_genotype_2$SNPname
      
      EAP_2 <- as.data.frame(correl_on_step(snp_mat = locus_genotype_val_2,
                                            expr_vec = as.numeric(as.character(unlist(expression_data_2[,-1]))),
                                            cortype = cortype))
      
      EAP_2$SNPname <- locus_genotype_2$SNPname
      EAP_2$ID <- locus_genotype_2$ID
      
      ###################################################
      
    }
    
    ##############_EQUALIZE_EAP1_AND_EAP2_#############
    
    gwgw_checked <- check_gwgw(DA1 = EAP_1, DA2 = EAP_2)
    gwgw_checked[[1]] <- ref_alt_check(locus_genotype = gwgw_checked[[2]],
                                       DA = gwgw_checked[[1]])
    
    ###################################################
    
  }
  
  # return(gwgw_checked)
  
  
  if (plot_type == "range"){
    
    plot_theta <- theta_range(gwgw_checked[[1]],gwgw_checked[[2]],
                              gene_1, gene_2,
                              tissue_1, tissue_2,
                              gene_positions,
                              draw_lines = draw_lines, two_axes = two_axes)
    
  } else {
    
    plot_theta <- theta_ass(gwgw_checked[[1]],gwgw_checked[[2]],
                            gene_1, gene_2,
                            tissue_1, tissue_2,
                            gene_positions, weight_cor)
    
    
  }
  
  return(plot_theta)
  
}

###################################################

###################_THETA_RANGE_###################

# https://stackoverflow.com/a/66055331

train_sec <- function(primary, secondary, na.rm = TRUE) {
  
  # Thanks Henry Holm for including the na.rm argument!
  from <- range(secondary, na.rm = na.rm)
  to   <- range(primary, na.rm = na.rm)
  # Forward transform for the data
  forward <- function(x) {
    rescale(x, from = from, to = to)
  }
  # Reverse transform for the secondary axis
  reverse <- function(x) {
    rescale(x, from = to, to = from)
  }
  list(fwd = forward, rev = reverse)
  
}

#### MAIN

theta_range <- function(EAP_1 = EAP_1, EAP_2 = EAP_2, 
                        gene_1 = gene_1, gene_2 = gene_2,
                        tissue_1 = tissue_1, tissue_2 = tissue_2,
                        gene_positions = gene_positions,
                        draw_lines = draw_lines,
                        two_axes = two_axes){
  
  ###################################################
  
  if (gene_1==gene_2){
    
    EAP_1$gene <- paste0(gene_1,"_tissue_1")
    EAP_2$gene <- paste0(gene_2,"_tissue_2")
    
  } else {
    
    EAP_1$gene <- gene_1
    EAP_2$gene <- gene_2
    
  }
  
  
  EAP_2 <- EAP_2[,colnames(EAP_1)]
  
  EAP_1$p_value <- (-log10(EAP_1$p_value))
  EAP_2$p_value <- (-log10(EAP_2$p_value))
  
  merged_data <- rbind(EAP_1, EAP_2)
  merged_data$gene <- as.factor(merged_data$gene)
  
  gene_positions$strand <- as.factor(gene_positions$strand)
  
  
  
  if (gene_1==gene_2){
    
    if (gsub("_tissue_1", "",levels(merged_data$gene)[1])%in%gene_positions$pid){
      
      pid_gene <- gene_positions[gene_positions$pid==gsub("_tissue_1", "",levels(merged_data$gene)[1]),
                                 "gid"]
      
      levels(merged_data$gene)[1] <- paste0(pid_gene, "_tissue_1")
      levels(merged_data$gene)[2] <- paste0(pid_gene, "_tissue_2")
      
    }
    
  } else {
    
    levels(merged_data$gene) <- sapply(levels(merged_data$gene), function(x){
      
      if (x%in%gene_positions$pid){
        
        x <- gene_positions[gene_positions$pid==x,"gid"]
        
      } else {x <- x}
      
      x
      
    })
    
    
  }
  
  
  if (gene_1%in%gene_positions$pid){
    
    gene_1_name <- gene_positions[gene_positions$pid==gene_1,"gid"]
    
  } else {gene_1_name <- gene_1}
  
  if (gene_2%in%gene_positions$pid){
    
    gene_2_name <- gene_positions[gene_positions$pid==gene_2,"gid"]
    
  } else {gene_2_name <- gene_2}
  
  merged_data$pos <- as.numeric(sapply(merged_data$ID, function(x){
    
    strsplit(x,":")[[1]][2]
    
  }))/1000
  
  
  
  myColors <- c("#f8766d","#00bfc4")
  names(myColors) <- levels(merged_data$gene)
  
  ###################################################
  
  if (isTRUE(two_axes)){
    
    ###################################################
    
    cmerged_data <- merged_data[1:nrow(EAP_1),]
    
    cmerged_data$p_value_sec <- EAP_2$p_value
    second_axis <- with(cmerged_data,train_sec(p_value,p_value_sec))
    
    plot_range <- ggplot(cmerged_data, aes(x = pos)) +
      geom_point(aes(y = p_value, 
                     colour = names(myColors)[names(myColors)==merged_data$gene[1]])) + 
      geom_point(aes(y = second_axis$fwd(p_value_sec), 
                     colour = names(myColors)[names(myColors)!=merged_data$gene[1]])) + 
      theme_bw() +
      scale_y_continuous(name = paste0("-log10(p-value), ", gene_1_name),
                         sec.axis = sec_axis(~second_axis$rev(.), 
                                             name = paste0("-log10(p-value), ", gene_2_name))) +
      labs(title = paste0("Colocolisation between AP of ", 
                          gene_1_name," in ", tissue_1,
                          " and ", gene_2_name," in ", tissue_2)) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 15)) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            plot.title = element_text(size=9)) +
      labs(colour = "Genes") +
      scale_color_manual(values = myColors, labels = names(myColors))
    
    ###################################################
    
  } else {
    
    ###################################################
    
    plot_range <- ggplot(merged_data, aes(x = pos, y = p_value, color = gene)) +
      geom_point() + theme_bw() +
      ylab("-log10(p-value)") +
      labs(title = paste0("Colocolisation between AP of ", 
                          gene_1_name," in ", tissue_1,
                          " and ", gene_2_name," in ", tissue_2)) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 15)) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            plot.title = element_text(size=9)) +
      scale_color_manual(values = myColors)
    
    ###################################################
    
  }
  
  ###################################################
  
  gen_on_plot <- as.data.frame(gene_positions[gene_positions$gid%in%gsub("_tissue_1",
                                                                         "",gsub("_tissue_2",
                                                                                 "",levels(merged_data$gene))),])
  
  gen_on_plot$start <- gen_on_plot$start/1000
  gen_on_plot$end <- gen_on_plot$end/1000
  
  for (i in 1:nrow(gen_on_plot)){
    
    if (!(min(merged_data$pos)<=gen_on_plot[i,"start"]&
          gen_on_plot[i,"start"]<=max(merged_data$pos))){
      
      gen_on_plot[i,c("start","end")] <- max(merged_data$pos)
      
    } else if (gen_on_plot[i,"end"]>max(merged_data$pos)) {
      
      gen_on_plot[i,"end"] <- max(merged_data$pos)
      
    }
    
  }
  
  ###################################################
  
  gen_on_plot_gg <- as.data.frame(list(c(gen_on_plot$start,gen_on_plot$end),
                                       c(gen_on_plot$gid,gen_on_plot$gid),
                                       c(gen_on_plot$strand,gen_on_plot$strand)))
  
  colnames(gen_on_plot_gg) <- c("pos","gid","strand")
  levels(gen_on_plot_gg$strand) <- levels(gene_positions$strand) 
  gen_on_plot_gg$gid <- as.factor(gen_on_plot_gg$gid)
  
  if (any(gen_on_plot$start!=gen_on_plot$end)){
    
    gen_on_plot_gg <- gen_on_plot_gg[gen_on_plot_gg$gid%in%
                                       gen_on_plot[gen_on_plot$start!=gen_on_plot$end,"gid"],]
    
    genes_plot <- ggplot(gen_on_plot_gg, aes(x = pos, y = strand)) +
      geom_point(aes(colour = gid), size = 3) + 
      geom_line(aes(group = gid, colour = gid), size = 4) +
      theme_bw() + geom_hline(yintercept = c("+","-"), 
                              color = "black", linetype = "dashed")+
      xlab("Position, kb") + ylab("Strand") +
      scale_y_discrete(limits =c("-","+")) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 15),
                         limits =c(min(merged_data$pos),max(merged_data$pos)),
                         labels = scales::number_format(accuracy = 1, decimal.mark = '.')) +
      theme(axis.text.x = element_text(angle = 15),
            axis.text.y = element_text(size = 18)) +
      scale_color_manual(values = myColors, drop = F)
    
  } else {
    
    genes_plot <- ggplot(gen_on_plot_gg, aes(x = pos, y = strand)) +
      theme_bw() + geom_hline(yintercept = c("+","-"), 
                              color = "black", linetype = "dashed")+
      geom_blank() +
      xlab("Position, kb") + ylab("Strand") +
      scale_y_discrete(limits =c("-","+")) +
      scale_x_continuous(breaks = scales::pretty_breaks(n = 15),
                         limits =c(min(merged_data$pos),max(merged_data$pos)),
                         labels = scales::number_format(accuracy = 1, decimal.mark = '.')) +
      theme(axis.text.x = element_text(angle = 15),
            axis.text.y = element_text(size = 18)) +
      scale_color_manual(values = myColors, drop = F)
    
  }
  
  
  
  if (!is.null(draw_lines)){
    
    plot_range <- plot_range + geom_vline(xintercept = draw_lines, 
                                          color = "blue", linetype = "dashed")
    
  }
  
  plot_theta <- ggarrange(plot_range, genes_plot, ncol=1, nrow=2, 
                          common.legend = TRUE, legend="bottom",
                          align = "v", heights = c(4,1))
  
  return(plot_theta)
  
  
}

###################################################

####################_THETA_ASS_####################

theta_ass <- function(EAP_1 = EAP_1, EAP_2 = EAP_2, 
                      gene_1 = gene_1, gene_2 = gene_2,
                      tissue_1 = tissue_1, tissue_2 = tissue_2,
                      gene_positions = gene_positions,
                      weight_cor = weight_cor){
  
  if (gene_1%in%gene_positions$pid){
    
    gene_1 <- gene_positions[gene_positions$pid==gene_1,"gid"]
    
  }
  
  if (gene_2%in%gene_positions$pid){
    
    gene_2 <- gene_positions[gene_positions$pid==gene_2,"gid"]
    
  }
  
  X_logp <- -log10(EAP_1$p_value)
  Y_logp <- -log10(EAP_2$p_value)
  
  X_beta <- EAP_1$cor_or_beta
  Y_beta <- EAP_2$cor_or_beta
  
  w_i <- mapply(max, X_logp/max(X_logp), Y_logp/max(Y_logp))
  
  X_w <- sum(X_logp*w_i)/sum(w_i)
  Y_w <- sum(Y_logp*w_i)/sum(w_i)
  
  X_sigma <- sqrt(sum(w_i*(X_logp-X_w)^2)/sum(w_i))
  Y_sigma <- sqrt(sum(w_i*(Y_logp-Y_w)^2)/sum(w_i))
  
  
  if (isTRUE(weight_cor)) {
    
    X_plot <- w_i*((X_logp-X_w)/X_sigma)
    Y_plot <- w_i*((Y_logp-Y_w)/Y_sigma)*sign(X_beta*Y_beta)
    
  } else {
    
    X_plot <- X_logp
    Y_plot <- Y_logp*sign(X_beta*Y_beta)
    
  }
  
  
  
  data_gg <- as.data.frame(cbind(X_plot,Y_plot))
  
  if (isTRUE(weight_cor)) {
    
    plot_theta <- ggplot(data_gg,aes(x = X_plot, y = Y_plot)) +
      geom_point() + theme_bw() +
      xlab(paste0("Weighted -log10(p-value), ", gene_1)) + 
      ylab(paste0("Weighted -log10(p-value), ", gene_2)) +
      labs(title = paste0("Correlation between weighted -log10(p-value) of ", 
                          gene_1," in ", tissue_1,
                          " and ", gene_2," in ", tissue_2)) +
      theme(plot.title = element_text(size=9))
    
  } else {
    
    plot_theta <- ggplot(data_gg,aes(x = X_plot, y = Y_plot)) +
      geom_point() + theme_bw() +
      xlab(paste0("-log10(p-value), ", gene_1)) + 
      ylab(paste0("-log10(p-value), ", gene_2)) +
      labs(title = paste0("Correlation between -log10(p-value) of ", 
                          gene_1," in ", tissue_1,
                          " and ", gene_2," in ", tissue_2)) +
      theme(plot.title = element_text(size=9))
    
  }
  
  
  
  return(plot_theta)
  
  
  
}

###################################################

####################_PRINT_HMM_####################

print_HMM_res <- function(GWAS_dat = GWAS_dat,
                          known = known,
                          predicted = predicted,
                          chr = chr){
  
  ## known - name of column with ground-true info or NULL
  ## predicted - name of column with HMM predicted labels or NULL
  
  GWAS_dat <- GWAS_dat[GWAS_dat$chr==chr,]
 
  
  if (!is.null(known)){
    
    loci_borders_known <- as.data.frame(matrix(ncol = 3,nrow = 0))
    colnames(loci_borders_known) <- c("chr","left","right")
    
    prev <- 0
    a <- 1
    
    for (i in 1:nrow(GWAS_dat)){
      
      if (GWAS_dat[[known]][i]==0&prev==0){} else
        if (GWAS_dat[[known]][i]==1&prev==0){
          
          prev <- 1
          loci_borders_known[a,1] <- as.character(GWAS_dat$chr[i])
          loci_borders_known[a,2] <- GWAS_dat$pos[i]
          
        } else if (GWAS_dat[[known]][i]==1&prev==1){} else
          if (GWAS_dat[[known]][i]==0&prev==1){
            
            prev <- 0
            loci_borders_known[a,3] <- GWAS_dat$pos[i-1]
            a <- a+1
            
          }
      
    }
    
  }
  
  if (!is.null(predicted)){
    
    loci_borders_predicted <- as.data.frame(matrix(ncol = 3,nrow = 0))
    colnames(loci_borders_predicted) <- c("chr","left","right")
    
    prev <- 0
    a <- 1
    
    for (i in 1:nrow(GWAS_dat)){
      
      if (GWAS_dat[[predicted]][i]==0&prev==0){} else
        if (GWAS_dat[[predicted]][i]==1&prev==0){
          
          prev <- 1
          loci_borders_predicted[a,1] <- as.character(GWAS_dat$chr[i])
          loci_borders_predicted[a,2] <- GWAS_dat$pos[i]
          
        } else if (GWAS_dat[[predicted]][i]==1&prev==1){} else
          if (GWAS_dat[[predicted]][i]==0&prev==1){
            
            prev <- 0
            loci_borders_predicted[a,3] <- GWAS_dat$pos[i-1]
            a <- a+1
            
          }
      
    }
    
  }
  
  GWAS_dat$fakepos <- c(1:nrow(GWAS_dat))
  
  gwas_peak <- ggplot(GWAS_dat, aes(x = fakepos, y = -log10(p_value))) +
    geom_point() + 
    #stat_binhex() +
    theme_bw() + 
    labs(title = "HMM predicted and known")+
    xlab("Position, kb") + ylab("-log10(p-value)") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 15)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 11)) +
    theme(axis.text.x = element_text(angle = 15, size=10),
          axis.text.y = element_text(size=12))
  
  if (!is.null(known)){
    
    loci_borders_known$fakeleft <- apply(loci_borders_known,1, function(x){
      
      GWAS_dat[GWAS_dat$chr==x[1]&GWAS_dat$pos==as.numeric(x[2]),"fakepos"]
      
    })
    
    loci_borders_known$fakeright <- apply(loci_borders_known,1, function(x){
      
      GWAS_dat[GWAS_dat$chr==x[1]&GWAS_dat$pos==as.numeric(x[3]),"fakepos",]
      
    })
    
    gwas_peak <- gwas_peak + geom_vline(xintercept = unlist(loci_borders_known[,4:5]), 
                                        color = "blue", linetype = "dashed")
    
  }
  
  if (!is.null(predicted)){
    
    loci_borders_predicted$fakeleft <- unlist(apply(loci_borders_predicted,1, function(x){
      
      GWAS_dat[GWAS_dat$chr==x[[1]]&GWAS_dat$pos==as.numeric(x[[2]]),"fakepos",][1]
      
    }))
    
    loci_borders_predicted$fakeright <- unlist(apply(loci_borders_predicted,1, function(x){
      
      GWAS_dat[GWAS_dat$chr==x[[1]]&GWAS_dat$pos==as.numeric(x[[3]]),"fakepos",][1]
      
    }))
    
    gwas_peak <- gwas_peak + geom_rect(data = loci_borders_predicted,
                                       mapping = aes(xmin = fakeleft, xmax = fakeright),
                                       ymin = 0, ymax = max(-log10(GWAS_dat$p_value)),
                                       fill = "red", alpha = 0.2, inherit.aes = F)
    
  }
  
  plot(gwas_peak)
  
  
}


###################################################

#####################_TESTRUNS_AREA_#################################




