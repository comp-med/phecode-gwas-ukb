########################################
## function to compute cis-pQTL 
## colocalisation 

coloc.cis.pQTL <- function(res, soma, gene, chr.s, g.e, g.s, ld, r.t=F){
  
  ## 'res'   -- data frame with regional GWAS results
  ## 'soma'  -- name of the protein target
  ## 'gene'  -- associated protein coding gene (sometimes one protein target maps to multiple genes)
  ## 'chr.s' -- chromosome of the protein coding gene
  ## 'g.s'   -- start position of the protein coding gene
  ## 'g.e'   -- end position of the protein coding gene
  ## 'ld'    -- ld matrix to track consistency of signals
  ## 'r.t'   -- whehter to report back merged summary statistics
  
  ## load package to do so
  require(coloc)
  
  ## rename INDELs (coded I/D in Fenland)
  res                            <- as.data.frame(res)
  res[, c("ALLELE0", "ALLELE1")] <- t(apply(res[, c("ALLELE0", "ALLELE1")], 1, function(k){
    if(nchar(k[1]) > nchar(k[2])){
      return(c("I", "D"))
    }else if(nchar(k[2]) > nchar(k[1])){
      return(c("D", "I"))
    }else{
      return(toupper(k))
    }
  }))
  
  ## create MarkerName to identify possible multi-allelic INDELs that are not well coded in the proteomics data  
  res$pQTL.MarkerName <- apply(res[, c("CHROM", "GENPOS", "ALLELE0", "ALLELE1")], 1, function(x){
      paste0("chr", as.numeric(x[1]), ":", as.numeric(x[2]), "_", paste(sort(x[3:4]), collapse = "_"))
  })
  
  ## discard ambiguous variants
  ii                  <- table(res$pQTL.MarkerName)
  res                 <- subset(res, pQTL.MarkerName %in% names(ii[ii == 1]))
  
  #-----------------------------------------#
  ##-- 	       import protein data       --##
  #-----------------------------------------#
  
  ## get the protein target - gene combination of interest
  res.soma    <- dir("<path to file>")
  res.soma    <- grep(soma, res.soma, value=T)
  res.soma    <- grep(gene, res.soma, value=T)
  res.soma    <- grep(g.s, res.soma, value=T)
  res.soma    <- grep(g.e, res.soma, value=T)
  
  ## read the relevant data (has been fine-mapped with SuSie before-hand)
  res.soma    <- fread(paste0("<path to file>", res.soma))
  
  ## identify top SNP by Z-score 
  top.snps    <- res.soma[!is.na(cs)]
  top.snps    <- top.snps[order(cs, -pip)]
  ## take only the strongest for each credible set
  top.snps[, ind := 1:.N, by="cs"]
  top.snps    <- top.snps[ind == 1]
  
  ## store the LD pattern across top SNPs (convert to data frame to ease downstream operations)
  ld.top.snps <- as.data.frame(res.soma)
  ld.top.snps <- lapply(1:nrow(top.snps), function(x){
    ## get all SNPs and corresponding LD
    tmp        <- paste0("R2.", x)
    tmp        <- ld.top.snps[which(ld.top.snps[, tmp] >= .8), c("MarkerName", "id", "rsid", tmp)]
    ## edit names
    names(tmp) <- c("MarkerName.proxy", "id.proxy", "rsid.proxy", "R2")
    ## add top SNP
    tmp        <- merge(as.data.frame(top.snps[x, c("MarkerName", "id", "rsid", "cs")]), tmp, suffix=c(".lead", ".proxy"))
    ## do some renaming to ease downstream coding
    names(tmp) <- c("MarkerName.lead", "id.lead", "rsid.lead", "cs", "MarkerName.proxy", "id.proxy", "rsid.proxy", "R2.proxy")
    ## return
    return(tmp)
  })
  ## combine everything
  ld.top.snps <- do.call(rbind, ld.top.snps)
  
  #-----------------------------------------#
  ##-- 	      combine both data sets     --##
  #-----------------------------------------#
  
  ## combine with results from SomaLogic
  res.all         <- merge(res, res.soma, by.x="pQTL.MarkerName", by.y="MarkerName", suffixes = c(".protein", ".trait"))
  ## convert to data frame
  res.all         <- as.data.frame(res.all)
  
  ## align effect estimates based on UKB results
  res.all$Effect  <- ifelse(toupper(res.all$Allele1) == res.all$ALLELE1, res.all$Effect, -res.all$Effect)
  ## delete everything no longer needed
  res.all$Allele1 <- res.all$Allele2 <- NULL
  
  #-----------------------------------------#
  ##-- 	          run coloc test         --##
  #-----------------------------------------#
  
  ## prepare input list (N.B.: use numeric SNP id to be able to map back to LD file)
  D1        <- list(beta=res.all$Effect, varbeta=res.all$StdErr^2, 
                    type="quant", 
                    N=max(res.all$TotalSampleSize, na.rm=T), 
                    sdY=1,
                    MAF=res.all$MAF,
                    snp=res.all$snp.id,
                    position=1:nrow(res.all))
  
  ## binary outcome
  D2          <- list(beta=res.all$BETA, varbeta=res.all$SE^2, 
                      type="cc",
                      s=.01, 
                      N=max(res.all$N),
                      MAF=res.all$MAF,
                      snp=res.all$snp.id,
                      position=1:nrow(res.all))
  
  ## run coloc 
  res.coloc   <- coloc.signals(D1, D2, method="single", p12=5e-6)
  
  #-----------------------------------------#
  ##--      cross-check lead signals     --##
  #-----------------------------------------#
  
  ## add LD with fine-mapped variants 
  res.coloc         <- as.data.frame(res.coloc$summary)
  
  ## compute ld between selected lead hits
  res.coloc$ld.top  <- apply(res.coloc[, c("hit1", "hit2"), drop=F], 1, function(k) ld[as.numeric(k[1]), as.numeric(k[2])]^2)
  
  #-----------------------------------------#
  ##--     add additional information    --##
  #-----------------------------------------#
  
  ## add effect estimates (choose SNP with highers posterior prob. for the shared signal, no regional sentinels taken forward!, only LD informatio)
  res.coloc            <- merge(res.coloc, res.all, by.x = "best4", by.y = "snp.id")
  
  ## add MarkerName for lead signals - protein
  res.coloc            <- merge(res.coloc, res.all[, c("pQTL.MarkerName", "snp.id")], by.x = "hit1", by.y = "snp.id", suffixes = c(".shared", ".protein.proxy"))
  ## same for phenotype
  res.coloc            <- merge(res.coloc, res.all[, c("pQTL.MarkerName", "snp.id")], by.x = "hit2", by.y = "snp.id", suffixes = c(".shared", ".phecode.proxy"))
  ## rename last edit
  names(res.coloc)     <- gsub("^pQTL.MarkerName$", "pQTL.MarkerName.phecode.proxy", names(res.coloc))
  
  ## add data on lead credible variant(s), may drop some without mapping lead credible variant
  res.coloc            <- merge(res.coloc, ld.top.snps[, c("MarkerName.proxy", "R2.proxy")], by.x="pQTL.MarkerName.protein.proxy", by.y="MarkerName.proxy", all.x=T, suffixes = c(".lead", ".overlap"))
  ## same for lead signal from the phecode GWAS
  res.coloc$ld.phecode <- ld[res.coloc$best2, res$snp.id[which.max(res[, "LOG10P"])]]^2
  
  ## add protein
  res.coloc$pheno      <- soma
  
  ## convert to data frame
  res.coloc            <- as.data.frame(res.coloc)
  
  #-----------------------------------------#
  ##--           return results          --##
  #-----------------------------------------#
  
  cat("\n\n----------------------\n")
  print(res.coloc)
  
  a <- c("pheno", "nsnps", "ld.top", "R2.proxy", "ld.phecode", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf", 
         "MarkerName", "rsid", "chr", "pos", "MAF", "ALLELE0", "ALLELE1",
         "Effect", "StdErr", "Pvalue", "TotalSampleSize", "BETA", "SE", "LOG10P", "N")
  
  print(a[!(a %in% names(res.coloc))])
  
  
  ## subset to what is really needed (N.B.: MarkerName is the one create with INDEL information from UK Biobank)
  res.coloc        <- res.coloc[, c("pheno", "nsnps", "ld.top", "R2.proxy", "ld.phecode", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf", 
                                    "MarkerName", "rsid", "chr", "pos", "MAF", "ALLELE0", "ALLELE1",
                                    "Effect", "StdErr", "Pvalue", "TotalSampleSize", "BETA", "SE", "LOG10P", "N")]
  
  ## change names to be more explanatory
  names(res.coloc) <- c("pheno", "nsnps", "ld.top.overlap", "ld.sens.protein", "ld.sens.phecode", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf",
                        "MarkerName", "rsid", "chr", "pos", "MAF", "ALLELE0", "ALLELE1",
                        "Effect.protein", "StdErr.protein", "Pvalue.protein", "N.protein", "Effect.trait", "StdErr.trait", "log10p.trait", "N.trait")
  
  if(r.t){
    return(list(res.coloc, res.all))
  }else{
    return(res.coloc)
  }
   
  
}