###########################################
## function to retrieve regional gene 
## dosage files to compute alleles

get.LD <- function(chr.s, pos.s, pos.e, pheno){
  
  ## 'chr.s' -- Chromosome
  ## 'pos.s' -- start position
  ## 'pos.e' -- end position
  ## 'pheno' -- include phecode to be able to process potentially overlapping intervals

  ## run as as a bash job
  print(paste("./scripts/get_dosages_LD_matrix.sh", chr.s, pos.s, pos.e, pheno))
  system(paste("./scripts/get_dosages_LD_matrix.sh", chr.s, pos.s, pos.e, pheno))
  
  ## read the dosage file
  require(data.table)
  tmp                 <- fread(paste("tmpdir/tmp", pheno, chr.s, pos.s, pos.e, "dosage", sep="."), sep=" ", header=T, data.table=F)
  ## transpose
  rownames(tmp)       <- tmp$SNPID
  ## store allele information to realign effect estimates afterwards
  tmp.info            <- tmp[, 1:6]
  tmp                 <- t(tmp[,-c(1:6)])
  ## re-transform to data set (keep X in mind for variable names)
  tmp                 <- data.frame(ID_1=rownames(tmp), tmp)
  ## make IDS numeric to ease merging afterwards
  tmp$ID_1            <- as.numeric(tmp$ID_1)
  
  ## create another column to info to map names from the SNP data set
  tmp.info$id         <- sapply(tmp.info$SNPID, function(x) ifelse(substr(x, 1, 2) == "rs", x, paste0("X", gsub(":", ".", x))))
  ## edit some IDs (X-chromosome)
  tmp.info$id         <- gsub("XX", "X", tmp.info$id)
  tmp.info$id         <- gsub("XAffx-", "Affx.", tmp.info$id)
  tmp.info$id         <- gsub(",", ".", tmp.info$id)
  ## set to those included in the data set, just to be sure
  tmp.info            <- subset(tmp.info, id %in% names(tmp))
  
  ## create MarkerName to merge to results files
  tmp.info$MarkerName <- apply(tmp.info, 1, function(x){
    paste0("chr", ifelse(x[1] == "X", 23, as.numeric(x[1])), ":", as.numeric(x[4]), "_", paste(sort(x[5:6]), collapse = "_"))
  })
  
  return(list(tmp, tmp.info))
  
}