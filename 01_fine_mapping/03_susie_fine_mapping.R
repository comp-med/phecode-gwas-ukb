#!/usr/bin/env Rscript

## script to perfrom fine-mapping for a given region in a GWAS
## Maik Pietzner
rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)
## avoid conversion of numbers
options(scipen = 1)
# print(R.Version())

setwd("<path to file>")

## --> import parameters <-- ##

## get regional coordinates
pheno <- args[1]
chr.s <- as.numeric(args[2])
pos.s <- as.numeric(args[3])
pos.e <- as.numeric(args[4])

## get regional coordinates
# tmp   <- read.table("input/phecode.fine.mapping.regions.large.txt")
# j     <- 1
# pheno <- tmp$V1[j]
# chr.s <- tmp$V2[j]
# pos.s <- tmp$V3[j]
# pos.e <- tmp$V4[j]

# pheno <- "date_702.1"
# chr.s <- 5
# pos.s <- 33644592
# pos.e <- 34251693


#-----------------------------------------#
##--     load regional assoc stats     --##
#-----------------------------------------#

## package for faster data handling
require(data.table)

## read the relevant data
res            <- paste0("zcat <path to file>/gwas_results/", pheno, ".allchr.results.gz | awk -v chr=", chr.s, " -v low=", pos.s, " -v upp=", pos.e, 
                         " '{if(($1 == chr && $2 >= low && $2 <= upp) || NR == 1) print $0}'")
res            <- data.table::fread(cmd = res, sep = " ", header = T, data.table = F)

## drop SNPs that have possibly failed in REGENIE
res            <- subset(res, is.na(EXTRA))

## create MarkerName to enable mapping to the LD file
res$MarkerName <- apply(res, 1, function(x){
  paste0("chr", as.numeric(x[1]), ":", as.numeric(x[2]), "_", paste(sort(x[4:5]), collapse = "_"))
})

#-----------------------------------------#
##--            LD-matrix              --##
#-----------------------------------------#

## write list of SNPs to be queried to file 
write.table(res$ID, paste("tmpdir/snplist", pheno, chr.s, pos.s, pos.e, "lst", sep="."), col.names = F, row.names = F, quote = F)

## snp data to be queried
tmp.z        <- res[, c("ID", "CHROM", "GENPOS", "ALLELE0", "ALLELE1")]
names(tmp.z) <- c("rsid", "chromosome", "position", "allele1", "allele2")

## adopt chromosome if needed
if(chr.s < 10){
  print("tada")
  tmp.z$chromosome <- paste0("0", tmp.z$chromosome)
}else if(chr.s == 23){
  tmp.z$chromosome <- "X"
}

## check for input
print(head(tmp.z))

## write to file
write.table(tmp.z, paste("tmpdir/snpz", pheno, chr.s, pos.s, pos.e, "z", sep="."), row.names = F, quote = F)

## --> create master file for LDstore2 <-- ##

## assign entries
m.file      <- data.frame(z=paste("tmpdir/snpz", pheno, chr.s, pos.s, pos.e, "z", sep="."),
                          bgen=paste("tmpdir/filtered", pheno, chr.s, pos.s, pos.e, "bgen", sep="."),
                          bgi=paste("tmpdir/filtered", pheno, chr.s, pos.s, pos.e, "bgen.bgi", sep="."),
                          ld=paste("tmpdir/ld", pheno, chr.s, pos.s, pos.e, "ld", sep="."),
                          incl="<path to file>",
                          n_samples=ifelse(chr.s != 23, 487409, 486757))

## write to file
write.table(m.file, paste("tmpdir/master", pheno, chr.s, pos.s, pos.e, "z", sep="."), sep=";", row.names = F, quote = F)

cat("--------------------------------------\n")
cat("computing LD matrix\n")

## submit the job
system(paste("./scripts/get_LD_matrix_ldstore.sh", chr.s, pos.s, pos.e, pheno))

## read in matrix
ld         <- fread(paste("tmpdir/ld", pheno, chr.s, pos.s, pos.e, "ld", sep="."), data.table = F)
## print for checking
print(ld[1:5,1:5])

## create identifier column in results to keep the mapping to the LD matrix
res$snp.id <- 1:nrow(res) 

cat("Done\n")
cat("--------------------------------------\n")

#-----------------------------------------#
##--          run fine-mapping         --##
#-----------------------------------------#

require('susieR')

## run fine mapping to obtain 95%-credible sets (throws an error due to possible rounding errors in LD matrix)
set.seed(42)
## run and catch possible errors
res.fine <- tryCatch(
  {
    susie_rss(res$BETA/res$SE, as.matrix(ld), L = 10, coverage = .95, min_abs_corr=.1, max_iter = 10000)
  }, error=function(e){
    return(list(pip=rep(NA, nrow(res)),
                sets=list(cs=NA),
                converged=F))
  })

## look at results
if(!is.na(res.fine$pip[1])){
  print(summary(res.fine))
}else{
  cat("-----------------------------\n")
  cat("SuSiE failed to find solution\n")
}

## do if condition in case now credible sets exist
if(length(res.fine$sets$cs) > 0 & res.fine$converged == T){
  
  #-----------------------------------------#
  ##--            store results          --##
  #-----------------------------------------#
  
  ## add PIP to results (may drop some rare SNPs with no LD information)
  tmp      <- data.frame(snp.id=gsub("V", "", names(res.fine$pip)), pip=res.fine$pip)
  res      <- merge(res, tmp, all.x = T, by="snp.id")
  
  ## add credible set information
  tmp      <- res.fine$sets$cs
  ## add info to the results file
  res$cs   <- NA
  
  #------------------------------------------------#
  ##-- run joint model to ensure GWAS threshold --##
  #------------------------------------------------#
  
  ## identify top SNP for each credible set
  top.snp         <- lapply(tmp, function(x){
    ## ensure that exact mapping is preserved
    foo <- res[which(res$snp.id %in% x), ]
    return(foo$snp.id[which.max(foo$pip)])
  })
  
  ## import label data set to identify phenotyp batch
  label           <- fread("<path to file>")
  ## identify the one of interest
  label           <- label[id == pheno]
  ## import phenotype
  phecode         <- fread(paste("<path to file>", label$batch, "txt", sep = "."),
                           select=c("IID", pheno))
  ## import covariates
  covs            <- fread("<path to file>")
  ## combine both
  phecode         <- merge(phecode, covs)
  ## import inclusion list of unrelated white European individuals
  inc.list        <- fread("<path to file>")
  ## subset to those
  phecode         <- phecode[IID %in% inc.list$ID_1]

  ## check for sex-specific outcomes
  if(label$sex == "Female"){
    phecode <- phecode[sex == 1]
  }
  if(label$sex == "Male"){
    phecode <- phecode[sex == 0]
  }
  
  ## import function to do so
  source("scripts/get_LD.R")
  
  ## write list of SNPs to be queried to file 
  write.table(res$ID[which(res$snp.id %in% unlist(top.snp))], paste("tmpdir/snplist", pheno, chr.s, pos.s, pos.e, "lst", sep="."), col.names = F, row.names = F, quote = F)
  
  ## import 
  foo              <- get.LD(chr.s, pos.s, pos.e, pheno)
  ## separate out into two data sets to ease downstream computation
  snp.dat          <- foo[[1]]
  snp.info         <- foo[[2]]
  ## delete and clear up
  rm(foo); gc()
  ## add results variant ID to snp.info
  snp.info         <- merge(snp.info, res[, c("ID", "snp.id", "MarkerName")], by.x="MarkerName", by.y="MarkerName")
  ## add credible sets
  # snp.info$cs      <- sapply(snp.info$snp.id, function(k) grep(k, tmp))
  snp.info$cs      <- sapply(snp.info$snp.id, function(k){
    ## go through each credible set
    ii <- lapply(tmp, function(z) k %in% z)
    return(which(ii == T))
  })
  
  ## add SNP data
  phecode          <- merge(phecode, snp.dat, by.x="IID", by.y="ID_1")
  
  #--------------------------------#
  ##--      run joint model     --##
  #--------------------------------#
  
  ## run model
  m.joint         <- summary(glm(paste(pheno, "~", paste(c(snp.info$id, "age", "sex", "centre", "batch", paste0("pc", 1:10)), collapse = " + ")), phecode, family = binomial))$coefficients[-1, ,drop=F]
  ## look at SNPs only
  m.joint         <- as.data.frame(m.joint[snp.info$id, ,drop=F])
  
  ## add LOG10P from the GWAS
  m.joint$id      <- row.names(m.joint) 
  m.joint         <- merge(m.joint, snp.info)
  m.joint         <- merge(m.joint, res[, c("snp.id", "LOG10P")])
  ## change one name
  names(m.joint)  <- gsub("Pr\\(>\\|z\\|\\)", "PVALUE.joint", names(m.joint))
  
  ## keep only what passes filter
  m.joint         <- subset(m.joint, LOG10P > 7.3 | PVALUE.joint < 5e-8)
  
  ## only proceed it al least one remaining
  if(nrow(m.joint) > 0){
    
    ## credible sets to be kept
    ii              <- m.joint$cs
    
    ## subset credible sets
    tmp             <- tmp[ii]
    
    ## very last round
    top.snp         <- lapply(tmp, function(x){
      ## ensure that exact mapping is preserved
      foo <- res[which(res$snp.id %in% x), ]
      return(foo$snp.id[which.max(foo$pip)])
    })

    ## subset to SNPs still included
    snp.info        <- subset(snp.info, snp.id %in% unlist(top.snp))
    ## re-assign credible set numbers
    snp.info$cs     <- sapply(snp.info$snp.id, function(k){
      ## go through each credible set
      ii <- lapply(tmp, function(z) k %in% z)
      return(which(ii == T))
    })
    
    ## loop through all credible sets
    for(j in 1:length(tmp)){
      
      ## set all variants to j if included in the current credible set
      ## careful needs proper redefinement of tmp, since it no longer contains the credible sets
      res$cs[which(res$snp.id %in% tmp[[j]])] <- j
      
      ## add LD with top variant in the set
      ii                                      <- res[which(res$cs == j),]
      ii                                      <- ii$snp.id[which.max(ii$pip)]
      ## add LD
      ii                                      <- data.frame(snp.id=1:nrow(ld), R2=ld[,ii]^2)
      ## rename
      names(ii)                               <- c("snp.id", paste0("R2.", j))
      res                                     <- merge(res, ii)
    }
    
    ## store the final credible sets
    write.table(res, paste("output/credible.set", pheno, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names = F)
    
    #--------------------------------#
    ##--      run joint model     --##
    #--------------------------------#
    
    ## run final joint model
    m.joint         <- summary(glm(paste(pheno, "~", paste(c(snp.info$id, "age", "sex", "centre", "batch", paste0("pc", 1:10)), collapse = " + ")), phecode, family = binomial))$coefficients[-1, ,drop=F]
    ## look at SNPs only
    m.joint         <- as.data.frame(m.joint[snp.info$id, , drop=F])
    ## keep what is needed
    names(m.joint) <- c("BETA.joint", "SE.joint", "TVAL.joint", "PVALUE.joint")
    ## add id
    m.joint$id     <- rownames(m.joint)
    ## add credible set identifier
    m.joint$cs     <- snp.info$cs
    ## add allele coding
    m.joint        <- merge(snp.info[, c("id", "MarkerName", "alleleA", "alleleB")], m.joint)
    ## combine with marginal stats
    m.joint        <- merge(res, m.joint)
    ## add method
    m.joint$method <- "SuSiE"
    
    ## store results
    write.table(m.joint, paste("joint_models//fine.mapped", pheno, chr.s, pos.s, pos.e, "txt", sep = "."), sep="\t", row.names=F)
    
  }else{
    
    cat("found no credible sets, use Wakefield approximation assuming only one causal variant\n")
    
    #-----------------------------------------#
    ##-- implement Wakefield approximation --##
    #-----------------------------------------#
    
    ## compute PIP probabilities (https://cran.r-project.org/web/packages/corrcoverage/vignettes/corrected-coverage.html)
    require(corrcoverage)
    z                        <- res$BETA/res$SE
    v                        <- res$SE^2
    res$pip                  <- ppfunc(z, v, W = .2)
    
    ## compute the 99%-credible set
    cred.set                 <- credset(res$pip, thr = .95)
    
    ## add to the data
    res$cs                   <- NA
    res$cs[cred.set$credset] <- 1
    
    ## add LD (this time only one causal variant)
    ii                       <- res$snp.id[which.max(res$pip)]
    ## add LD
    ii                       <- data.frame(snp.id=1:nrow(ld), R2.1=ld[,ii]^2)
    ## rename
    res                      <- merge(res, ii)
    
    ## store the final credible sets
    write.table(res, paste("output/credible.set", pheno, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names = F)
    
    #------------------------------------------#
    ##-- pseudo joint model to ease merging --##
    #------------------------------------------#
    
    ## write list of SNPs to be queried to file 
    write.table(res$ID[which.max(res$pip)], paste("tmpdir/snplist", pheno, chr.s, pos.s, pos.e, "lst", sep="."), col.names = F, row.names = F, quote = F)
    
    ## import 
    foo              <- get.LD(chr.s, pos.s, pos.e, pheno)
    ## separate out into two data sets to ease downstream computation
    snp.dat          <- foo[[1]]
    snp.info         <- foo[[2]]
    ## delete and clear up
    rm(foo); gc()
    ## add results variant ID to snp.info
    snp.info         <- merge(snp.info, res[, c("ID", "snp.id", "MarkerName")], by.x="MarkerName", by.y="MarkerName")
    ## add credible sets
    snp.info$cs      <- 1
    
    ## delete prevoius SNPs from the data set
    phecode          <- phecode[, c("IID", pheno, "FID", "age", "sex", "centre", "batch", paste0("pc", 1:10)), with=F]
    ## add SNP data
    phecode          <- merge(phecode, snp.dat, by.x="IID", by.y="ID_1")
    
    ## --> run joint model <-- ##
    
    ## run final joint model
    m.joint         <- summary(glm(paste(pheno, "~", paste(c(snp.info$id, "age", "sex", "centre", "batch", paste0("pc", 1:10)), collapse = " + ")), phecode, family = binomial))$coefficients[-1, ,drop=F]
    ## look at SNPs only
    m.joint         <- as.data.frame(m.joint[snp.info$id, , drop=F])
    ## keep what is needed
    names(m.joint) <- c("BETA.joint", "SE.joint", "TVAL.joint", "PVALUE.joint")
    ## add id
    m.joint$id     <- rownames(m.joint)
    ## add credible set identifier
    m.joint$cs     <- snp.info$cs
    ## add allele coding
    m.joint        <- merge(snp.info[, c("id", "MarkerName", "alleleA", "alleleB")], m.joint)
    ## combine with marginal stats
    m.joint        <- merge(res, m.joint)
    ## add method
    m.joint$method <- "Wakefield"
    
    ## store results
    write.table(m.joint, paste("joint_models//fine.mapped", pheno, chr.s, pos.s, pos.e, "txt", sep = "."), sep="\t", row.names=F)
    
  }
  
}else{
  
  cat("found no credible sets, use Wakefield approximation assuming only one causal variant\n")
  
  #-----------------------------------------#
  ##-- implement Wakefield approximation --##
  #-----------------------------------------#
  
  ## compute PIP probabilities (https://cran.r-project.org/web/packages/corrcoverage/vignettes/corrected-coverage.html)
  require(corrcoverage)
  z                        <- res$BETA/res$SE
  v                        <- res$SE^2
  res$pip                  <- ppfunc(z, v, W = .2)
  
  ## compute the 99%-credible set
  cred.set                 <- credset(res$pip, thr = .95)
  
  ## add to the data
  res$cs                   <- NA
  res$cs[cred.set$credset] <- 1
  
  ## add LD (this time only one causal variant)
  ii                       <- res$snp.id[which.max(res$pip)]
  ## add LD
  ii                       <- data.frame(snp.id=1:nrow(ld), R2.1=ld[,ii]^2)
  ## rename
  res                      <- merge(res, ii)
  
  ## store the final credible sets
  write.table(res, paste("output/credible.set", pheno, chr.s, pos.s, pos.e, "txt", sep="."), sep="\t", row.names = F)
  
  #------------------------------------------#
  ##-- pseudo joint model to ease merging --##
  #------------------------------------------#
  
  ## import label data set to identify phenotyp batch
  label           <- fread("<path to file>")
  ## identify the one of interest
  label           <- label[id == pheno]
  ## import phenotype
  phecode         <- fread(paste("<path to file>", label$batch, "txt", sep = "."),
                           select=c("IID", pheno))
  ## import covariates
  covs            <- fread("<path to file>")
  ## combine both
  phecode         <- merge(phecode, covs)
  ## import inclusion list of unrelated white European individuals
  inc.list        <- fread("<path to file>")
  ## subset to those
  phecode         <- phecode[IID %in% inc.list$ID_1]
  
  ## check for sex-specific outcomes
  if(label$sex == "Female"){
    phecode <- phecode[sex == 1]
  }
  if(label$sex == "Male"){
    phecode <- phecode[sex == 0]
  }
  
  ## import function to do so
  source("scripts/get_LD.R")
  
  ## write list of SNPs to be queried to file 
  write.table(res$ID[which.max(res$pip)], paste("tmpdir/snplist", pheno, chr.s, pos.s, pos.e, "lst", sep="."), col.names = F, row.names = F, quote = F)
  
  ## import 
  foo              <- get.LD(chr.s, pos.s, pos.e, pheno)
  ## separate out into two data sets to ease downstream computation
  snp.dat          <- foo[[1]]
  snp.info         <- foo[[2]]
  ## delete and clear up
  rm(foo); gc()
  ## add results variant ID to snp.info
  snp.info         <- merge(snp.info, res[, c("ID", "snp.id")], by.x="rsid", by.y="ID")
  ## add credible sets
  snp.info$cs      <- 1
  
  ## add SNP data
  phecode          <- merge(phecode, snp.dat, by.x="IID", by.y="ID_1")
  
  ## --> run joint model <-- ##
  
  ## run final joint model
  m.joint         <- summary(glm(paste(pheno, "~", paste(c(snp.info$id, "age", "sex", "centre", "batch", paste0("pc", 1:10)), collapse = " + ")), phecode, family = binomial))$coefficients[-1, ,drop=F]
  ## look at SNPs only
  m.joint         <- as.data.frame(m.joint[snp.info$id, , drop=F])
  ## keep what is needed
  names(m.joint) <- c("BETA.joint", "SE.joint", "TVAL.joint", "PVALUE.joint")
  ## add id
  m.joint$id     <- rownames(m.joint)
  ## add credible set identifier
  m.joint$cs     <- snp.info$cs
  ## add allele coding
  m.joint        <- merge(snp.info[, c("id", "MarkerName", "alleleA", "alleleB")], m.joint)
  ## combine with marginal stats
  m.joint        <- merge(res, m.joint)
  ## add method
  m.joint$method <- "Wakefield"
  
  ## store results
  write.table(m.joint, paste("joint_models//fine.mapped", pheno, chr.s, pos.s, pos.e, "txt", sep = "."), sep="\t", row.names=F)
}

## delete LD file
system(paste("rm tmpdir/*", pheno, chr.s, pos.s, pos.e, "*", sep="."))

cat("-----------------------\n")
cat("DONE!")