######################################
## function to run coloc across GTEx
## tissues

coloc.eQTL <- function(res, chr.s, pos.s, pos.e, gene, ld, r.t=F){
  
  ## 'res'     -- summary stats for the outcome in pGWAS format style
  ## 'chr.s'   -- chromosome
  ## 'pos.s'   -- start position
  ## 'pos.e'   -- end position
  ## 'gene'    -- ENSEMBL identfier for the gene of interest
  ## 'ld'      -- LD-matrix 
  ## 'r.t'     -- whether or not to return summary stats
  
  #-----------------------------------------#
  ##--        import tissue data         --##
  #-----------------------------------------#
  
  ## get all GTEx files
  ii          <- dir("<path to file>/GTEx_v8_Science_2020/eQTL") 
  
  ## get unqiue tissues
  tissue      <- grep("csv", ii, value=T)
  tissue      <- unique(sapply(tissue, function(x) strsplit(x, "\\.")[[1]][1]))
  
  ## do in parallel
  require(doMC)
  registerDoMC(6)

  ## get the Stats
  res.gtex    <- mclapply(tissue, function(x){
    ## get all information for the respective gene
    tmp        <- paste0("grep ", gene, " ", "<path to file>/GTEx_v8_Science_2020/eQTL/",
                         x, ".v8.EUR.allpairs.chr", ifelse(chr.s == 23, "X", chr.s),".csv")
    # print(tmp)
    tmp        <- data.table::fread(cmd=tmp, sep=",", header=F, data.table = F)
    ## proceed only if gene is available in the respective tissue
    if(nrow(tmp)>1){
      ## add header
      names(tmp) <- c("nr", "phenotype_id", "variant_id", "tss_distance", "maf", "ma_samples", "ma_count", "pval_nominal", "slope", "slope_se")
      ## add tissue
      tmp$tissue <- x
      ## compute N
      tmp$N      <- (tmp$ma_count/tmp$maf)/2
      ## keep only finite estiamtes
      tmp        <- subset(tmp, is.finite(slope))
      return(tmp[, c("phenotype_id", "variant_id", "N", "slope", "slope_se", "pval_nominal", "tissue")])
    }else{
      cat(gene, "not found in", x, "\n")
    }
    
  }, mc.cores=6)
  
  ## --> proceed only if found in at least one tissue <-- ##
  if(!is.null(res.gtex[[1]])){
    
    ## combine into one data frame
    res.gtex <- do.call(rbind, res.gtex)
    
    print(dim(res.gtex))
    
    ## reshape
    library(tidyr)
    res.gtex <- pivot_wider(res.gtex, id_cols = c("variant_id", "phenotype_id"),
                            values_from = c("slope", "slope_se", "pval_nominal", "N"),
                            names_from = "tissue", names_sep = ".")
    ## convert to data frame
    res.gtex <- data.frame(res.gtex)
    
    ## subset to tissues with available gene expression
    tissue   <- grep("slope_se", names(res.gtex), value=T)
    tissue   <- gsub("slope_se\\.", "", tissue)
    
    ## create position and allele columns
    res.gtex$pos.hg38 <- as.numeric(sapply(res.gtex$variant_id, function(x) strsplit(x, "_")[[1]][2]))
    ## coding allele
    res.gtex$Allele1  <- sapply(res.gtex$variant_id, function(x) strsplit(x, "_")[[1]][4])
    res.gtex$Allele2  <- sapply(res.gtex$variant_id, function(x) strsplit(x, "_")[[1]][3])
    
    ## delete non-biallelic variants
    ii                <- table(res.gtex$pos.hg38)
    ii                <- names(ii[which(ii > 1)])
    if(length(ii) > 0){
      res.gtex <- res.gtex[-which(res.gtex$pos %in% ii),]
    }
    
    #-----------------------------------------#
    ##--         combine everything        --##
    #-----------------------------------------#
    
    ## combine
    res.all                <- merge(res.gtex, res, by="pos.hg38", suffixes = c(".gtex", ".pheno"))
    
    ## proceed only if enough observations+
    if(nrow(res.all) >= 300){

      ## create aligned effect estimate for phenotype
      res.all$Effect.aligned <- ifelse(res.all$Allele1 == res.all$ALLELE1, res.all$BETA, -res.all$BETA)
      
      ## omit NAs
      res.all$EXTRA          <- NULL
      res.all                <- na.omit(res.all)
      
      print(dim(res.all))
      
      ## sort
      res.all                <- res.all[order(res.all$pos.hg38),]
      
      ## compute minor allele frequency
      res.all$MAF            <- ifelse(res.all$A1FREQ >= .5, 1-res.all$A1FREQ, res.all$A1FREQ)
      
      ## store the top hit for the protein
      top.snp                <- which.max(res$LOG10P)
      top.snp                <- res$snp.id[top.snp]
      
      ## get new top snp in the merged data set
      top.merged             <- which.max(res.all$LOG10P)
      top.merged             <- res.all$snp.id[top.merged]
      
      ## compute the LD between both
      ld.sens                <- ld[top.merged, top.snp]^2
      
      #-----------------------------------------#
      ##--            run coloc              --##
      #-----------------------------------------#
      
      require(coloc)
      
      ## run coloc in parallele
      res.coloc <- mclapply(tissue, function(x){
        
        print(x)
        
        #-----------------------------------------#
        ##-- 	         sanity check            --##
        #-----------------------------------------#
        
        ## get the top SNP for the outcome
        io     <- res.all$snp.id[which.max(abs(res.all[, paste0("slope.",x)]/res.all[, paste0("slope_se.",x)]))]
        
        ## keep names
        ld.top <- ld[io, top.merged]^2
        
        #-----------------------------------------#
        ##-- 	            run coloc            --##
        #-----------------------------------------#
        
        ## prepare input for protein
        D1          <- list(beta=res.all$Effect.aligned, 
                            varbeta=res.all$SE^2, 
                            type="cc", 
                            N=max(res.all$N), 
                            s=.01,
                            MAF=res.all$MAF,
                            snp=res.all$snp.id,
                            position=1:nrow(res.all))
        
        ## prepare input for GTEx
        D2          <- list(beta=res.all[, paste0("slope.",x)], 
                            varbeta=res.all[, paste0("slope_se.",x)]^2,
                            # pvalues=res.all[, paste0("pval_nominal.",x)],
                            type="quant",
                            N=max(res.all[, paste0("N.", x)]),
                            MAF=res.all$MAF,
                            snp=res.all$snp.id,
                            sdY=1,
                            position=1:nrow(res.all))
        
        ## do naive coloc as well
        naive.coloc <- coloc.signals(D1, D2, method="single", p12=5e-6)
        
        ## add checks to the data
        naive.coloc$summary$ld.check.top  <- ld.top
        naive.coloc$summary$ld.check.sens <- ld.sens
        
        ## add the trait id
        naive.coloc$summary$tissue        <- x
        
        ## write results to file
        naive.coloc                       <- naive.coloc$summary
        # print(naive.coloc)
        
        ## give back data for top pQTL
        ii                                <- which(res.all$snp.id == top.merged) 
        
        ## add top SNP from phenotype with estimates
        naive.coloc                       <- data.frame(naive.coloc, 
                                                        res.all[ii, c("phenotype_id", "CHROM", "pos.hg38", "GENPOS", "variant_id", "MarkerName", "ID",
                                                                      "Allele1", "Allele2", "MAF", "Effect.aligned", "SE", "LOG10P",
                                                                      paste(c("slope", "slope_se", "pval_nominal"), x, sep="."))])
        ## edit names
        names(naive.coloc)                <- gsub(paste0("\\.", x),  "", names(naive.coloc))
        return(naive.coloc)
        
      }, mc.cores = 6)
      res.coloc <- do.call(rbind, res.coloc)
      
      if(r.t){
        return(list(res.coloc, res.all))
      }else{
        return(res.coloc)
      }
      
    }else{
      return(NULL)
    }
    
        
  }else{
    return(NULL)
  }
  
  
}