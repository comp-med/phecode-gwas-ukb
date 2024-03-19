######################################################
#### Perform cis-pQTL colocalisation              ####
#### Maik Pietzner                                ####
######################################################

rm(list=ls())
setwd("<path to file>")
options(stringsAsFactors = F)
load(".RData")

## packages needed
require(data.table)

###############################################
####      obtain fine-mapped results       ####
###############################################

## import fine-mapped results
res.fine <- fread("<path to file>")

## drop MHC region
res.fine <- res.fine[group != "6_11"]

## create new regional boundaries for each genetic variant
res.fine[, region_pQTL_end := GENPOS + 5e5] 
res.fine[, region_pQTL_start := pmax(0,GENPOS - 5e5)] 

###############################################
####          perpare for coloc            ####
###############################################

## export file with input parameters for coloc
write.table(res.fine[, c("id", "CHROM", "GENPOS", "region_pQTL_start", "region_pQTL_end")], "input.pQTL.pipeline.txt", sep="\t", row.names=F, col.names=F, quote=F)

###############################################
####           collate results             ####
###############################################

## gather output from joint models
ii               <- dir("../output/")

## set options before
options(scipen = 1)

#----------------------------------#
##--  check for possible empty  --##
#----------------------------------#

## import file size
f.size <- fread(cmd="ls -l ../output", header=F)
tail(sort(f.size$V5, decreasing = T), 200)
## only very few files failed due to missing gene mappings

#----------------------------------#
##--        import results      --##
#----------------------------------#

## faster reading
require(doMC)
registerDoMC(6)

## loop through all instances and collate basic information on credible sets
res.proteins <- mclapply(ii, function(x){
  
  ## import the correct stats file
  tmp                   <- fread(paste0("../output/", x), header=T, sep="\t", data.table=F)
  
  ## some might not have produces results due to missing overlapping genes
  if(nrow(tmp) > 0){
    
    ## do some annotations to be able to map back to olink genes
    x                     <- strsplit(x, "\\.")[[1]]
    ## add pheno, region start and end
    tmp$id                <- paste(unique(c(x[2], x[length(x)-4])), collapse = ".")
    tmp$region_pQTL_start <- as.numeric(x[length(x)-2])
    tmp$region_pQTL_end   <- as.numeric(x[length(x)-1])
    
    ## return
    return(tmp)
  }
  
}, mc.cores = 6)
## delete NULL entries
res.proteins <- plyr::compact(res.proteins)
## combine
res.proteins <- do.call(plyr::rbind.fill, res.proteins)
## do some QC
res.proteins <- subset(res.proteins, ld.sens.phecode > .8)
## subset to medium to high confidence findings
res.proteins <- subset(res.proteins, PP.H4.abf > .5)
## convert to data table to speed up queries
res.proteins <- as.data.table(res.proteins)

###############################################
####  create summary column for each locus ####
###############################################

## go through each fine-mapping result and add putative gene annotations
foo <- mclapply(1:nrow(res.fine), function(x){
  
  ## get the gene sets of most interest
  tmp <- res.proteins[ id == res.fine$id[x] & chr == res.fine$CHROM[x] & region_pQTL_start == res.fine$region_pQTL_start[x] & region_pQTL_end == res.fine$region_pQTL_end[x]]
  
  ## proceed only if any
  if(nrow(tmp) > 0){
    
    ## summarize eQTL results (by gene)
    tmp <- lapply(unique(tmp$ensembl_gene_id), function(k){
      ## subset
      foo <- tmp[ensembl_gene_id == k]
      ## return summary
      return(paste0(foo$hgnc_symbol[1], "(", sprintf("%.1f", max(foo$PP.H4.abf*100)), "%)-", nrow(foo), "-", k, "-", foo$Target[1]))
    })
    ## return results
    return(data.frame(res.fine[x,], pQTL_summary=paste(tmp, collapse = "||")))
  }else{
    return(res.fine[x,])
  }
  
}, mc.cores=6)
## combine everything
foo <- do.call(plyr::rbind.fill, foo)
## how many with at least suggestive evidence
nrow(subset(foo, !is.na(pQTL_summary)))
## prioritized 291 phecode - loci associations...

## rename
res.coloc <- foo

## write to file
write.table(res.coloc, "<path to file>", sep="\t", row.names=F)