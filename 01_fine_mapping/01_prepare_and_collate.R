######################################################
#### Fine-mapping of phecode GWAS results         ####
  #### Maik Pietzner                              ####
######################################################

rm(list=ls())
setwd("<path to file>")
options(stringsAsFactors = F)
load(".RData")

## --> packages needed <-- ##
require(data.table)

#######################################
####    prepare files for input    ####
#######################################

## header for results files
hd          <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA", "region_start", "region_end")

## collate results files from regional clumping
res.regions <- dir("<path to file>/regional_results/")
res.regions <- lapply(res.regions, function(x){
  
  ## import the file of interest
  if(file.size(paste0("<path to file>/regional_results/", x)) > 0){
    
    ## read the results file
    tmp         <- read.table(paste0("<path to file>/regional_results/", x), sep=" ", header = F)
    ## add names
    names(tmp)  <- hd
    ## add phenotype
    tmp$id      <- gsub("_regional_sentinels.txt", "", x)
    ## careful, replace TRUE for T in allele coding
    tmp$ALLELE0 <- gsub("TRUE", "T", tmp$ALLELE0)
    tmp$ALLELE1 <- gsub("TRUE", "T", tmp$ALLELE1)
    
    ## collapse potentially adjacent regions
    tmp         <- tmp[order(tmp$CHROM, tmp$region_start), ]
    
    ## return annotated results
    return(tmp)
  }else{
    cat("no results for", x, "\n")
  }
  
})
res.regions <- do.call(rbind, res.regions)
res.regions <- as.data.table(res.regions)
res.regions <- res.regions[order(CHROM, region_start)]

#-------------------------------------------#
##-- enlarge regions across all phecodes --##
#-------------------------------------------#

## create dummy boundaries to enforce at least some overlap 
## for regions across diseases
res.regions$region_start_cap <- res.regions$region_start + 2e5
res.regions$region_end_cap   <- res.regions$region_end - 2e5
## sort again
res.regions                  <- res.regions[order(CHROM, region_start_cap)]

## from https://stackoverflow.com/questions/28938147/how-to-flatten-merge-overlapping-time-periods
require(dplyr)
tmp      <- res.regions %>%
                ## sort
                arrange(CHROM, region_start_cap) %>% 
                group_by(CHROM) %>%
                ## create index for regions
                mutate(indx = c(0, cumsum(as.numeric(lead(region_start_cap)) >
                                            cummax(as.numeric(region_end_cap)))[-n()])) %>%
                group_by(CHROM, indx) %>%
                summarise(region_start_cap = first(region_start_cap), region_end_cap = last(region_end_cap)) %>% as.data.table()
## N = 1104

## set keys to indicate on what to define new regions on
setkey(tmp, CHROM, region_start_cap, region_end_cap)
## add information to regional results
foo                  <- foverlaps(res.regions, tmp, by.x=c("CHROM", "region_start_cap", "region_end_cap"))
## the prefix 'i' denotes the former regional boundaries before collapsing

## create new index for regions
foo$group            <- paste(foo$CHROM, foo$indx, sep="_")

## expand regions, that may have missed the signal of interest (close to the beginning of the chromosome)
foo[, ind := ifelse(GENPOS >= region_start_cap & GENPOS <= region_end_cap, 1, 0)]
foo[, region_start_cap := ifelse(GENPOS <= region_start_cap, GENPOS - 25e4, region_start_cap)]
foo[, region_end_cap := ifelse(GENPOS >= region_end_cap, GENPOS + 25e4, region_end_cap)]
## test again for not covered signals
foo[, ind := ifelse(GENPOS >= region_start_cap & GENPOS <= region_end_cap, 1, 0)]
foo[ind == 0]
## edit negative margins
foo[, region_start_cap := ifelse(region_start_cap < 0, 0, region_start_cap)]

## expand regions again, to cover at least 1MB
foo                  <- lapply(unique(foo$group), function(x){
  
  ## get all hits found in the respective region
  tmp                       <- foo[group == x]
  
  ## get new min and max boundaries
  tmp$region_start_collated <- min(tmp$region_start_cap)
  tmp$region_end_collated   <- max(tmp$region_end_cap)
  
  ## return
  return(tmp)
  
  
})
## combine again
foo                <- do.call(rbind, foo)
## width of each region
foo$width          <- foo$region_end_collated - foo$region_start_collated
## drop what is no longer needed
foo$region_end_cap <- foo$region_start_cap <- foo$i.region_end_cap <- foo$i.region_start_cap <- NULL

## re-assign
res.regions        <- foo

#-------------------------------------------#
##--        create input for SuSiE       --##
#-------------------------------------------#

## write to file 
tmp  <- unique(res.regions[group != "6_11", c("id", "CHROM", "region_start_collated", "region_end_collated", "width")])

## split jobs by width
tmp1 <- tmp[width < 3e6]
tmp2 <- tmp[width >= 3e6]
## write to file
options(scipen = 1)
write.table(tmp1, "phecode.fine.mapping.regions.small.txt", sep="\t", col.names = F, row.names = F, quote = F)
write.table(tmp2, "phecode.fine.mapping.regions.large.txt", sep="\t", col.names = F, row.names = F, quote = F)

#######################################
####         import results        ####
#######################################

## gather output from joint models
ii               <- dir("../joint_models//")

## set options before
options(scipen = 1)

#----------------------------------#
##--        import results      --##
#----------------------------------#

## faster reading
require(doMC)
registerDoMC(6)

## loop through all instances and collate basic information on credible sets
res.fine <- mclapply(ii, function(x){
  
  ## import the correct stats file
  tmp                       <- fread(paste0("../joint_models/", x), header=T, sep="\t")
  ## do some annotations to be able to map back to olink genes
  x                         <- strsplit(x, "\\.")[[1]]
  ## add pheno, region start and end
  tmp$id                    <- paste(unique(c(x[3], x[length(x)-4])), collapse = ".")
  tmp$region_start_collated <- as.numeric(x[length(x)-2])
  tmp$region_end_collated   <- as.numeric(x[length(x)-1])
  
  ## return
  return(tmp)
  
}, mc.cores = 6)
## combine
res.fine <- do.call(plyr::rbind.fill, res.fine)
## convert to data table
res.fine <- as.data.table(res.fine)

#---------------------------------------#
##-- combine with regional sentinels --##
#---------------------------------------#

## select strongest signal for each newly defined region (will collapse some stats)
res.regions <- res.regions[order(id, CHROM, region_start_collated)]
res.regions[, ind_region := 1:.N, by=c("id", "CHROM", "region_start_collated")]

## take only the strongest forward
foo         <- res.regions[ind_region == 1]

## create identifier
foo[, locus_id := paste(id, CHROM, region_start_collated, region_end_collated, sep = "|")]

## create identifier
res.fine[, locus_id := paste(id, CHROM, region_start_collated, region_end_collated, sep = "|")]

## --> check whether the lead signal from regional clumping is preserved <-- ##

## needs digging for LD proxies
require(doMC)
registerDoMC(6)

## get all fine-mapping results
ii <- dir("../output/")

## loop through all loci
foo <- mclapply(1:nrow(foo), function(x){
  ## variant ID
  v.id <- foo$ID[x]
  ## locus ID
  l.id <- foo$locus_id[x]
  ## import fine-mapping results
  res  <- paste("credible.set", gsub("\\|", ".", l.id), "txt", sep=".")
  ## check whether included in any of the credible sets
  if(res %in% ii){
    ## import fine-mapping results
    res  <- fread(paste0("../output/", res))
    ## subset to credible sets
    res   <- res[!is.na(cs)]
    v.inc <- ifelse(v.id %in% res$ID, T, F)
  }else{
    v.inc <- F
  }
  ## test whether the exact same variant is included
  return(data.table(foo[x,], included.fine.mapping = v.inc))
}, mc.cores=6)
## combine again
foo <- do.call(rbind, foo)
## look how many differed
table(foo$included.fine.mapping)
## N = 311 are discordant 
table(foo[ group != "6_11"]$included.fine.mapping)
## N = 79 outside of the MHC region, keep those for now, assuming that they decompose the marginal signal

## --> combine regional stats (MHC) with fine-mapping results <-- ##

## create MarkerName
foo$MarkerName <- apply(foo, 1, function(x){
  paste0("chr", as.numeric(x[1]), ":", as.numeric(x[3]), "_", paste(sort(x[5:6]), collapse = "_"))
})

## create combined data set
tmp1           <- as.data.frame(res.fine)
tmp2           <- as.data.frame(foo[group == "6_11", intersect(names(res.fine), names(foo)), with=F])
res.all        <- plyr::rbind.fill(tmp1, tmp2)
## add locus group
res.all        <- merge(res.all, foo[, c("locus_id", "group")])
## convert to data table again
res.all        <- as.data.table(res.all)

#--------------------------#
##--      add label     --##
#--------------------------#

## import label
label          <- fread("<path to file>/Label.phecode.GWAS.20220301.txt")
res.all        <- merge(label[, c("id", "phecode", "phenotype", "sex", "category", "batch", "y.pos")], res.all)

#######################################
####      store for processing     ####
#######################################

## write to file
write.table(res.all, "Results.phecode.GWAS.fine.mapping.20220426.txt", sep="\t", row.names=F)

###############################################
####  create summary of pleiotrop regions  ####
###############################################

## get all regions
pleio.regions <- data.frame(table(res.all$group))
## add coordinates
pleio.regions <- merge(unique(res.all[, c("group", "CHROM", "region_start_collated", "region_end_collated")]), pleio.regions, 
                       by.x="group", by.y="Var1")
## write to file
write.table(pleio.regions, "Phecode.GWAS.region.count.20220624.txt", sep="\t", row.names=F)
