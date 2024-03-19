######################################################
#### Annotate variants using OpenTargets          ####
#### Maik Pietzner                                ####
######################################################

rm(list=ls())
setwd("<path to file>")
options(stringsAsFactors = F)
## avoid problems with scientific notation
options(scipen = 1)
load(".RData")

## packages needed
require(data.table)

###############################################
#### obtain list of variants to be queried ####
###############################################

## import fine-mapping results
res.fine <- fread("<path to file>")

###########################################################
####  create mapping table to allow query for build 38 ####
###########################################################

#-------------------------------------#
##--        create mapping         --##
#-------------------------------------#

## get all unique variants from regional results
res.var          <- unique(res.fine[, c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1")])
res.var          <- as.data.table(res.var)
## add unique A1FREQ (differs depending on whether selected for sex-specific outcomes or not)
res.var$A1FREQ   <- sapply(res.var$ID, function(x) res.fine$A1FREQ[which(res.fine$ID == x)][1]) 

## import 
foo              <- fread("<path to file>", sep="\t", header=T)

## create MarkerName based on build37 to merge 
foo[, MarkerName:=paste0("chr", as.numeric(gsub("X", 23, chr_id_b37)), ":", as.numeric(position_b37), "_", pmin(alt_allele,ref_allele), "_", pmax(alt_allele,ref_allele))]
res.var[, MarkerName:=paste0("chr", as.numeric(gsub("X", 23, CHROM)), ":", as.numeric(GENPOS), "_", pmin(ALLELE0,ALLELE1), "_", pmax(ALLELE0,ALLELE1))]

## add to variants of interest (increase first and drop allele mismatches afterwards)
## careful: some variants didn't match
res.var          <- merge(res.var, foo, by="MarkerName", all.x=T)
## non-unique variant mappings

## drop possible duplications by MAF mapping (take the one closer to MAF in UKBB)
res.var$maf.ukbb <- ifelse(res.var$A1FREQ > .5, 1-res.var$A1FREQ, res.var$A1FREQ)
res.var$maf.open <- ifelse(res.var$gnomad_nfe > .5, 1-res.var$gnomad_nfe, res.var$gnomad_nfe)
## replace NA for variants not in OpenTargets data base
res.var$maf.open[which(is.na(res.var$maf.open))] <- 0
## create difference
res.var$maf.diff <- abs(res.var$maf.ukbb - res.var$maf.open)
res.var          <- res.var[order(MarkerName, maf.diff)]
## create indice
res.var[, ind:=1:.N, by="MarkerName"]
res.var          <- res.var[ind == 1]
## create eventual ID to overlay with OT data
res.var$open.id  <- paste(res.var$chr_id, res.var$position, res.var$ref_allele, res.var$alt_allele,  sep="_")

#---------------------------------------#
##-- find proxies, if no exact match --##
#---------------------------------------#

## how many variants with no mathcing OT variant
tmp              <- res.var[is.na(alt_allele)]

## create comprehensive proxy matrix
require(doMC)
registerDoMC(6)

## get all credible set results
ii               <- dir("<path to file>")

## go through each locus and collect all proxies for lead signals (exclude MHC region)
proxy.list       <- mclapply(unique(res.fine[group != "6_11"]$locus_id), function(x){
  
  ## set option
  
  ## import fine-mapping results
  res         <- fread(paste("<path to file>", gsub("\\|", ".", x), "txt", sep="."))
  
  ## get all variants of interest
  top.snps    <- res[!is.na(cs)]
  top.snps    <- top.snps[order(cs, -pip)]
  ## take only the strongest for each credible set
  top.snps[, ind := 1:.N, by="cs"]
  top.snps    <- top.snps[ind == 1]
  top.snps    <- as.data.frame(top.snps)
  
  # store the LD pattern across top SNPs (convert to data frame to ease downstream operations)
  ld.top.snps <- as.data.frame(res)
  ld.top.snps <- lapply(1:nrow(top.snps), function(k){
    ## get all SNPs and corresponding LD
    tmp        <- paste0("R2.", k)
    tmp        <- ld.top.snps[which(ld.top.snps[, tmp] >= .1), c("MarkerName", "ID", tmp), drop=F]
    ## edit names
    names(tmp) <- c("MarkerName.proxy", "ID.proxy", "R2")
    print(tmp)
    ## add top SNP
    tmp        <- merge(top.snps[k, c("MarkerName", "ID"), drop=F], tmp, suffix=c(".lead", ".proxy"))
    ## do some renaming to ease downstream coding
    names(tmp) <- c("MarkerName.lead", "ID.lead", "MarkerName.proxy", "ID.proxy", "R2.proxy")
    ## return
    return(tmp)
  })
  ## combine everything
  ld.top.snps <- do.call(rbind, ld.top.snps)

  ## return
  return(ld.top.snps)
  
}, mc.cores=6)
## combine
proxy.list              <- do.call(rbind, proxy.list)
## drop possible duplications
proxy.list              <- unique(proxy.list)
## simplify
proxy.list$MarkerName.1 <- pmin(proxy.list$MarkerName.lead, proxy.list$MarkerName.proxy) 
proxy.list$MarkerName.2 <- pmax(proxy.list$MarkerName.lead, proxy.list$MarkerName.proxy) 
## drop what is no longer needed
proxy.list              <- unique(proxy.list[, c("MarkerName.1", "MarkerName.2", "R2.proxy")])

## store
write.table(proxy.list, "<path to file>", sep="\t", row.names=F)

## --> create sub-list of variants needed to be queried <-- ##

## subset to what would be minimally needed
tmp.proxy               <- subset(proxy.list, MarkerName.1 %in% tmp$MarkerName | MarkerName.2 %in% tmp$MarkerName)

## import SNP list to get allele frequency as well
tmp.snps                <- fread("zcat <path to file>/date_401.1.allchr.results.gz")
## generate MarkerName
tmp.snps[, MarkerName:=paste0("chr", as.numeric(gsub("X", 23, CHROM)), ":", as.numeric(GENPOS), "_", pmin(ALLELE0,ALLELE1), "_", pmax(ALLELE0,ALLELE1))]
## subset
tmp.snps                <- tmp.snps[MarkerName %in% tmp.proxy$MarkerName.1 | MarkerName %in% tmp.proxy$MarkerName.2]
## reduce to what is really needed
tmp.snps                <- tmp.snps[, c("MarkerName", "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ"), with=F]
## add OT annotation
tmp.snps                <- merge(tmp.snps, foo, by="MarkerName", all.x=T)
## drop non-matches
tmp.snps                <- tmp.snps[ !is.na(chr_id_b37)]
## drop possible duplication by MAF mapping (take the one closer to MAF in UKBB)
tmp.snps$maf.ukbb       <- ifelse(tmp.snps$A1FREQ > .5, 1-tmp.snps$A1FREQ, tmp.snps$A1FREQ)
tmp.snps$maf.open       <- ifelse(tmp.snps$gnomad_nfe > .5, 1-tmp.snps$gnomad_nfe, tmp.snps$gnomad_nfe)
## create difference
tmp.snps$maf.diff       <- abs(tmp.snps$maf.ukbb - tmp.snps$maf.open)
tmp.snps                <- tmp.snps[order(MarkerName, maf.diff)]
## create indice
tmp.snps[, ind:=1:.N, by="MarkerName"]
tmp.snps                <- tmp.snps[ind == 1]

## create eventual ID to overlay with OT data
tmp.snps$open.id        <- paste(tmp.snps$chr_id, tmp.snps$position, tmp.snps$ref_allele, tmp.snps$alt_allele,  sep="_")

###########################################################
####               query gene assignment               ####
###########################################################

## load packages to query
library(ghql)
library(jsonlite)
library(doMC)

## register multiple kernels to be used
registerDoMC(6)

## go through all entries
tmp <- mclapply(1:nrow(res.var), function(x){
  
  ## get the variant to be queried
  snp <- res.var$open.id[x]
  
  print(snp)
  
  # create GraphQL client
  cli <- GraphqlClient$new(
    url = "https://api.genetics.opentargets.org/graphql"
  )
  
  # create query object
  qry <- Query$new()
  qry$query('my_query', 'query getGenesForVariant($varId: String!){
                                genesForVariant(variantId: $varId) {
                                  variant, overallScore, 
                                  gene {
                                    id, symbol
                                       }
                                      }
                                    }')
  
  
  ## get the results
  tmp     <- fromJSON(cli$exec(qry$queries$my_query, list(varId=snp)), flatten = TRUE)$data$genesForVariant
  
  if(!is.null(nrow(tmp))){
    ## drop zeros
    tmp     <- subset(tmp, overallScore > 0)
    ## sort by overall score
    tmp     <- tmp[order(-tmp$overallScore),]
    ## create summary column
    tmp$scl <- apply(tmp[, c("gene.symbol", "gene.id", "overallScore")], 1, function(k){
      paste(k[1], k[2], sprintf("%.4f", as.numeric(k[3])), sep="_")
    }) 
    
    ## return information needed
    return(data.frame(res.var[x,], v2g.top=tmp$scl[1], v2g.top.score=tmp$overallScore[1], v2g.all=paste(tmp$scl, collapse = "|")))
    
  }else{
    
    # return information needed
    return(data.frame(res.var[x,], v2g.top="", v2g.top.score=NA, v2g.all=""))
  }
  

}, mc.cores = 6)
tmp     <- do.call(rbind, tmp)
## re-assign
res.v2g <- tmp

#----------------------------------------#
##--      same for proxy variants     --##
#----------------------------------------#

## go through all entries
tmp <- mclapply(1:nrow(tmp.snps), function(x){
  
  ## get the variant to be queried
  snp <- tmp.snps$open.id[x]
  
  print(snp)
  
  # create GraphQL client
  cli <- GraphqlClient$new(
    url = "https://api.genetics.opentargets.org/graphql"
  )
  
  # create query object
  qry <- Query$new()
  qry$query('my_query', 'query getGenesForVariant($varId: String!){
                                genesForVariant(variantId: $varId) {
                                  variant, overallScore, 
                                  gene {
                                    id, symbol
                                       }
                                      }
                                    }')
  
  
  ## get the results
  tmp     <- fromJSON(cli$exec(qry$queries$my_query, list(varId=snp)), flatten = TRUE)$data$genesForVariant
  
  if(!is.null(nrow(tmp))){
    ## drop zeros
    tmp     <- subset(tmp, overallScore > 0)
    ## sort by overall score
    tmp     <- tmp[order(-tmp$overallScore),]
    ## create summary column
    tmp$scl <- apply(tmp[, c("gene.symbol", "gene.id", "overallScore")], 1, function(k){
      paste(k[1], k[2], sprintf("%.4f", as.numeric(k[3])), sep="_")
    }) 
    
    ## return information needed
    return(data.frame(tmp.snps[x,], v2g.top=tmp$scl[1], v2g.top.score=tmp$overallScore[1], v2g.all=paste(tmp$scl, collapse = "|")))
    
  }else{
    
    # return information needed
    return(data.frame(tmp.snps[x,], v2g.top="", v2g.top.score=NA, v2g.all=""))
  }
  
  
}, mc.cores = 6)
tmp     <- do.call(rbind, tmp)
## drop what didn't map
tmp     <- subset(tmp, !is.na(v2g.top.score))

#----------------------------------------#
##--      Fuse wherever possible      --##
#----------------------------------------#

## do in parallel
registerDoMC(6)

## define header to be kept
ii <- c("MarkerName", "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ")
## define smaller head to be added
jj <- c("open.id", "v2g.top", "v2g.top.score", "v2g.all", "gene_id_any", "gene_id_any_distance", "gene_id_prot_coding", "gene_id_prot_coding_distance", "most_severe_consequence", "phred")

## go through and do the assignment
foo <- mclapply(1:nrow(res.v2g), function(x){
  
  ## get the current variant of interest
  res <- res.v2g[x,]
  
  ## test whether assignment is missing
  if(is.na(res$chr_id_b37)){
    ## get all possible proxies of the selected variant
    v.proxy <- subset(tmp.proxy, MarkerName.1 == res$MarkerName | MarkerName.2 == res$MarkerName)
    ## proceed only if any
    if(nrow(v.proxy) > 0){
      ## look which of those map
      v.proxy      <- subset(v.proxy, MarkerName.1 %in% tmp$MarkerName | MarkerName.2 %in% tmp$MarkerName)
      ## proceed only if any mapping
      if(nrow(v.proxy) > 0){
        ## take only the strongest one forward
        v.proxy    <- v.proxy[which.max(v.proxy$R2.proxy), c("MarkerName.1", "MarkerName.2")]
        ## get the required stats
        v.proxy    <- v.proxy[v.proxy != res$MarkerName]
        v.proxy    <- tmp[which(tmp$MarkerName == v.proxy),]
        ## return information
        return(data.table(res.v2g[x, ii], v.proxy[, jj], by.proxy=v.proxy$MarkerName))
      }else{
        ## return what has already been mapped
        return(data.table(res.v2g[x, c(ii,jj)], by.proxy="none_found"))
      }
    }else{
      ## return what has already been mapped
      return(data.table(res.v2g[x, c(ii,jj)], by.proxy="none_found"))
    }
  }else{
    ## return what has already been mapped
    return(data.table(res.v2g[x, c(ii,jj)], by.proxy="no"))
  }
}, mc.cores=6)
## combine everything
foo     <- do.call(rbind, foo)

## --> save as final list of assignment <-- ##
res.v2g <- foo

###########################################################
####          query phenotypic associations            ####
###########################################################

## get all json files from v2d
ii <- dir("../tmpdir/")
ii <- grep("json", ii, value=T)

## do in parallel
registerDoMC(6)

## loop over all files
res.v2d <- mclapply(ii, function(x){
  
  ## import file
  tmp          <- ndjson::stream_in(paste0("../tmpdir/", x))
  ## create identifier to be able to merge with results
  tmp[, open.id := paste(tag_chrom, tag_pos, tag_ref, tag_alt, sep="_")]
  ## keep only what has been found in the GWAS (makes already use of proxies identified in the previous step)
  tmp          <- subset(tmp, open.id %in% res.v2g$open.id)
  ## do some filtering based on R2 and missing values
  tmp          <- tmp[!is.na(overall_r2) | lead_pos == tag_pos]
  ## apply LD filtering 
  tmp          <- tmp[overall_r2 >= .5 | lead_pos == tag_pos]
  ## drop possible duplication
  tmp          <- tmp[order(open.id, study_id, -overall_r2)]
  tmp[, ind:=1:.N, by=c("open.id", "study_id")]
  ## keep only most
  tmp          <- tmp[ind == 1]
  ## keep only what is of immediate interest
  tmp          <- tmp[, c("ancestry_initial.0", "lead_chrom", "lead_pos", "lead_ref", "lead_alt", "open.id", "overall_r2", "trait_reported", "study_id",
                          "beta", "beta_ci_lower", "beta_ci_upper", "pval", "odds_ratio", "pval_exponent", "pval_mantissa", "source", "trait_category", "pmid")]
  gc()
  return(tmp)

}, mc.cores = 6)
res.v2d <- do.call(rbind, res.v2d)

## combine with info from our GWAS
foo     <- merge(res.fine, res.v2g[, -which(names(res.var) == "A1FREQ"), with=F])
foo     <- merge(foo, res.v2d, by="open.id", all.x=T, allow.cartesian = T)
res.v2d <- as.data.table(foo)

###########################################################
####    create collated list for fine-mapped stats     ####
###########################################################

## do in parallel
registerDoMC(12)

## collate information from different tables
res.fine.ot <- mclapply(1:nrow(res.fine), function(x){
  
  ## get common identifier to match to other tables
  snp <- res.fine$MarkerName[x]
  
  ## get the relevant V2G data
  v2g <- subset(res.v2g, MarkerName == snp)
  
  ## compress the relevant V2D data
  v2d <- res.v2d[MarkerName == snp]
  
  ## compress
  v2d <- data.frame(open.id=paste(sort(unique(v2d$open.id)), collapse = "||"),
                    trait_reported=paste(sort(unique(v2d$trait_reported)), collapse = "||"),
                    study_id=paste(sort(unique(v2d$study_id)), collapse = "||"),
                    source_ot=paste(sort(unique(v2d$source)), collapse = "||"),
                    num_reported=sum(!is.na(v2d$trait_reported)))
  
  ## report back 
  return(data.frame(res.fine[x,], v2g[, c("v2g.top", "v2g.top.score", "v2g.all")], v2d))
  
})
res.fine.ot <- do.call(rbind, res.fine.ot)
## convert to data table
res.fine.ot <- as.data.table(res.fine.ot)

## add column whether the variant isn't included in OpenTargets
res.fine.ot$ot.missing <- ifelse(res.fine.ot$open.id == "NA_NA_NA_NA", "yes", "no")

## how many, that haven't been seen yet
nrow(res.fine.ot[source_ot == "" & ot.missing == "no"])
View(res.fine.ot[source_ot == "" & ot.missing == "no"])

## write to file
write.table(res.fine.ot, "<path to file>", sep="\t", row.names=F)