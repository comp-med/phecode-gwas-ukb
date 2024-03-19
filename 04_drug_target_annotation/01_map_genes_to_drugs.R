######################################################
#### map candidate genes to drug targets using OT ####
#### Maik Pietzner                                ####
######################################################

rm(list=ls())
setwd("<path to file>")
options(stringsAsFactors = F)
load(".RData")

## packages needed
require(data.table)
require(doMC)
require(igraph)
require(rentrez)
require(httr)
require(jsonlite)

###########################################
####    import fine-mapping results    ####
###########################################

## import fine-mapping and annotation results from the phecode GWAS project
res.all  <- fread("<path to file>")
## import phecode information, including information on possibly matching drugs
res.ldsc <- fread("<path to file>")

###########################################
####      import drug information      ####
###########################################

## get all files in the directory
ii <- dir("<path to file>/download/molecule/")
ii <- grep("part", ii, value=T)

## run across multiple cores
registerDoMC(10)

## run in parallel
ot.drugs <- mclapply(ii, function(x){
  
  ## import each json data set and flatten
  tmp <- ndjson::stream_in(paste0("<path to file>/download/molecule/", x))
  
}, mc.cores = 10)
## combine all into one large data frame
ot.drugs <- do.call(plyr::rbind.fill, ot.drugs)
## convert to data table
ot.drugs <- as.data.table(ot.drugs) 

#-----------------------------#
##--   map OT ID to drugs  --##
#-----------------------------#

## create common disease column
ot.drugs$linkedDiseases.collated <- apply(ot.drugs[, grep("linkedDiseases.rows", names(ot.drugs), value=T), with=F], 1, function(x){
  ## compress
  paste(sort(unique(na.omit(x))), collapse = "|")
}) 

## create common target column
ot.drugs$linkedTargets.collated <- apply(ot.drugs[, grep("linkedTargets.rows", names(ot.drugs), value=T), with=F], 1, function(x){
  ## compress
  paste(sort(unique(na.omit(x))), collapse = "|")
}) 

###########################################
####  convert candidate causal genes   ####
###########################################

## use BiomaRt
require(biomaRt)

## create a list of gene symbols to be queried
genes                                   <- unique(unlist(lapply(res.all$candidate.gene.locus.r2, function(x) strsplit(x, "\\|")[[1]])))

## get data on build 37
gene.ensembl                            <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37) 

## obtain a list of all protein coding genes
tmp.genes                               <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id','external_gene_name', 'gene_biotype'),
                                                 filters = c('external_gene_name'),
                                                 values = list(genes),
                                                 mart = gene.ensembl)

## add entry to GWAS results
res.all$candidate.gene.locus.r2.ensembl <- sapply(res.all$candidate.gene.locus.r2, function(x){
  
  ## get all genes
  tmp <- strsplit(x, "\\|")[[1]]
  ## report all matching ensembl entries back
  return(paste(sort(unique(tmp.genes$ensembl_gene_id[which(tmp.genes$external_gene_name %in% tmp)])), collapse = "|"))
  
})

###################################
####  import diseases from OT  #### 
###################################

## get all files in the directory
ii <- dir("<path to file>/download/diseases/")
ii <- grep("part", ii, value=T)

## run across multiple cores
registerDoMC(10)

## run in parallel
ot.diseases <- mclapply(ii, function(x){
  
  ## import each json data set and flatten
  tmp <- ndjson::stream_in(paste0("<path to file>/download/diseases/", x))
  
}, mc.cores = 10)
## combine all into one large data frame
ot.diseases <- do.call(plyr::rbind.fill, ot.diseases)
## convert to data table
ot.diseases <- as.data.table(ot.diseases)
## drop soma names
ot.diseases <- ot.diseases[, grep(c("descendants|children"), names(ot.diseases), value = T, invert = T), with=F]

## create compressed columns for IDs
ot.diseases$dbXRefs.collated <- apply(ot.diseases[, grep("dbXRefs", names(ot.diseases), value=T), with=F], 1, function(x){
  ## compress
  paste(sort(unique(na.omit(x))), collapse = "|")
}) 
## omit dots from ICD10 codes
ot.diseases$dbXRefs.collated <- gsub("\\.", "", ot.diseases$dbXRefs.collated)

###################################
####    import MoA of drugs    #### 
###################################

## get all files in the directory
ii <- dir("<path to file>/download/mechanismOfAction/")
ii <- grep("part", ii, value=T)

## package to read parquet files
require(arrow)

## run across multiple cores
registerDoMC(10)

## run in parallel
ot.drugs.moa <- mclapply(ii, function(x){
  
  ## import each json data set and flatten
  tmp <- read_parquet(paste0("<path to file>/download/mechanismOfAction/", x), 
                      as_data_frame = T)
  ## drop what is not needed
  tmp <- tmp[, -which(names(tmp) == "references")]
  ## find list entries
  jj  <- unlist(lapply(tmp, is.list))
  ## go through and collapse lists
  for(j in names(jj[jj])){
    tmp[,j] <- unlist(lapply(tmp[[j]], function(k) paste0(k, collapse = "|")))
  }
  return(as.data.frame(tmp))
}, mc.cores = 10)
## combine into one large file
ot.drugs.moa <- do.call(rbind, ot.drugs.moa)

## unfold chemblIDs to ensure valid matching
ot.drugs.moa <- lapply(1:nrow(ot.drugs.moa), function(x){
  ## unfold
  id <- strsplit(ot.drugs.moa$chemblIds[x], "\\|")[[1]]
  ## return
  return(data.table(id=id, ot.drugs.moa[x,]))
})
## combine
ot.drugs.moa <- do.call(rbind, ot.drugs.moa)

###########################################
####        map drugs to variants      ####
###########################################

## map for each candidate genes potential drug targets
tmp <- mclapply(1:nrow(res.all), function(x){
  
  ## ids to search for
  gene  <- res.all$candidate.gene.locus.r2.ensembl[x]
  
  ## only if matching ENSEMBL ids
  if(gene != ""){
    ## search for ids in the OT drug target database
    drugs <- ot.drugs[grep(gene, ot.drugs$linkedTargets.collated),]
    
    ## return information
    return(data.table(res.all[x,],
                      drug.id=paste(drugs$id, collapse = "|"), 
                      drug.num=nrow(drugs),
                      drug.phase=paste(drugs$maximumClinicalTrialPhase, collapse = "|")))
  }else{
    ## return information
    return(data.table(res.all[x,], drug.id="", drug.num=0, drug.phase=""))
  }
  
}, mc.cores=10)
## combine again
tmp <- do.call(rbind, tmp)

## add information on whether there already is a treatment in trials
tmp <- merge(tmp, res.ldsc[, c("id", "drug.id")], by="id", suffixes = c(".gwas", ".trial"))

################################################
####          explore drug matchings        ####
################################################

## reassign
res.all                          <- tmp

## how often was the 'GWAS' drug already evaluated for the indication of interest
res.all$drug.gwas.trial.matching <- apply(res.all[, c("drug.id.gwas", "drug.id.trial"), with=T], 1, function(x){
  
  ## get all IDs
  gwas  <- strsplit(x[1], "\\|")[[1]]
  trial <- strsplit(x[2], "\\|")[[1]]
  
  ## report intersection
  return(paste(intersect(gwas, trial), collapse = "|"))
  
})

#----------------------------------#
##-- some numbers for the draft --##
#----------------------------------#

## create indicator by gene to avoid double counting
res.all <- res.all[order(id, candidate.gene.locus.r2, -candidate.gene.score, -LOG10P)]
res.all[, id.gene.ind := 1:.N, by=c("id", "candidate.gene.locus.r2")]

## how many combinations with at least one matching drug
nrow(res.all[ id.gene.ind == 1])
nrow(res.all[ drug.num > 0 & id.gene.ind == 1]) ## N = 887 out of 5683

## how many unique gene - drug assignments
length(unique(res.all[ drug.num > 0 & id.gene.ind == 1, candidate.gene.locus.r2])) ## N = 409
## how many genes in total
length(unique(res.all[ id.gene.ind == 1, candidate.gene.locus.r2])) ## N = 3088
## how many druggable genes
length(unique(unlist(lapply(ot.drugs$linkedTargets.collated, function(x) strsplit(x, "\\|")[[1]])))) ## N = 1462
## how many genes in the genome --> 19,817  Genome Res. 22, 1760–1774 (2012)

## test enrichment
(409/3088)/(1462/19817)
fisher.test(matrix(c(409,3088,1462,19817),2,2))
fisher.test(matrix(c(409,3088,1462,19817),2,2))$p.value

#-----------------------------------------#
##-- create table with full assignment --##
#-----------------------------------------#

## create new table
drug.table <- res.all[id.gene.ind == 1 & drug.num > 0, c("id", "phenotype", "locus_id", "gwas", "unreported", "candidate.gene.locus.r2", "candidate.gene.locus.r2.ensembl",  
                                                         "candidate.gene.score", "drug.id.gwas", "drug.num", "drug.id.trial", "drug.gwas.trial.matching")]
## extend to have each drug as separate entry
drug.table <- lapply(1:nrow(drug.table), function(x){
  
  ## get entry of interest
  tmp                  <- drug.table[x, , drop=F]
  ## report back all possible combinations
  tmp                  <- data.table(tmp, drug.id=strsplit(tmp$drug.id.gwas, "\\|")[[1]])
  ## add whether the drug is matching the indication
  jj                   <- strsplit(tmp$drug.gwas.trial.matching, "\\|")[[1]]
  tmp$match.indication <- tmp$drug.id %in% jj
  ## add number of gwas loci matching
  tmp$num.gwas.loci    <- nrow(res.all[ id == tmp$id[1] & candidate.gene.locus.r2 == tmp$candidate.gene.locus.r2[1]])
  return(tmp)
  
})
## combine
drug.table <- do.call(rbind, drug.table)

## add more information on drugs from OT
drug.table <- merge(drug.table, ot.drugs[, c("id", "blackBoxWarning", "description", "drugType", "hasBeenWithdrawn", "isApproved", "name", "synonyms", "tradeNames", "maximumClinicalTrialPhase", "linkedDiseases.collated", "linkedTargets.collated")],
                    by.x = "drug.id", by.y = "id")

################################################
####  augment drug matching using networks  ####
################################################

## expand gene annotations using 'signor 3.0' (https://signor.uniroma2.it/) Tue Sep 13 15:31:57 2022
sig.net  <- fread("Jul2022_release_signor_3_0.txt")

## obtain mapping of genes to uniprot IDs
tmp.prot <- unique(c(sig.net$IDA, sig.net$IDB))
## get mapping to gene names
tmp.prot <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id','external_gene_name', 'gene_biotype', 'uniprotswissprot'),
                  filters = c('uniprotswissprot'),
                  values = list(tmp.prot),
                  mart = gene.ensembl)

## add this information to sig.net
sig.net  <- merge(sig.net, tmp.prot, by.x = "IDA", by.y = "uniprotswissprot", all.x=T) ## careful: induces duplications
sig.net  <- merge(sig.net, tmp.prot, by.x = "IDB", by.y = "uniprotswissprot", all.x=T, suffixes = c(".A", ".B")) ## and nother round of duplications
sig.net  <- as.data.table(sig.net)

## add another column to res.all that augments candidate causal genes
tmp      <- lapply(1:nrow(res.all), function(x){
  
  ## get the relevant candidate causal genes
  genes <- strsplit(res.all$candidate.gene.locus.r2[x], "\\|")[[1]] 
  
  ## get the subset of merging entries in sig.net
  tmp   <- sig.net[ external_gene_name.A %in% genes | external_gene_name.B %in% genes ]
  
  ## report back
  return(data.table(res.all[x,], 
                    candidate.gene.locus.r2.net=paste(sort(unique(c(tmp$external_gene_name.A, tmp$external_gene_name.B))), collapse = "|"),
                    candidate.gene.locus.r2.net.ensembl=paste(sort(unique(c(tmp$ensembl_gene_id.A, tmp$ensembl_gene_id.B))), collapse = "|")))
  
})
## combine everything again
tmp <- do.call(rbind, tmp)

#-----------------------------------#
##--     rerun OT assignment     --##
#-----------------------------------#

## map for each candidate genes potential drug targets
tmp <- mclapply(1:nrow(tmp), function(x){
  
  ## ids to search for
  gene  <- tmp$candidate.gene.locus.r2.net.ensembl[x]
  
  ## only if matching ENSEMBL ids
  if(gene != ""){
    ## search for ids in the OT drug target database
    drugs <- ot.drugs[grep(gene, ot.drugs$linkedTargets.collated),]
    
    ## return information
    return(data.table(tmp[x,],
                      drug.id.net=paste(drugs$id, collapse = "|"), 
                      drug.num.net=nrow(drugs)))
  }else{
    ## return information
    return(data.table(tmp[x,], drug.id.net="", drug.num.net=0))
  }
  
}, mc.cores=10)
## combine again
tmp                              <- do.call(rbind, tmp)

## define intersection with possible targets already trialed for this condition
tmp$drug.gwas.net.trial.matching <- apply(tmp[, c("drug.id.net", "drug.id.trial"), with=T], 1, function(x){
  
  ## get all IDs
  gwas  <- strsplit(x[1], "\\|")[[1]]
  trial <- strsplit(x[2], "\\|")[[1]]
  
  ## report intersection
  return(paste(intersect(gwas, trial), collapse = "|"))
  
})

## how many with matching entries
nrow(tmp[ drug.gwas.net.trial.matching != "" & id.gene.ind == 1]) ## n = 387 with GWAS evidence
nrow(tmp[ drug.gwas.net.trial.matching != "" & id.gene.ind == 1 & candidate.gene.score >= 2]) ## n = 221 with GWAS evidence

## re-assign data table
res.all <- as.data.table(tmp)

#-----------------------------------------#
##-- create table with full assignment --##
#-----------------------------------------#

## create new table
drug.table.net <- res.all[id.gene.ind == 1 & (drug.num > 0 | drug.num.net > 0), c("id", "phenotype", "locus_id", "locus_id.r2", "gwas", "unreported", "candidate.gene.locus.r2", "candidate.gene.locus.r2.ensembl",  
                                                                                  "candidate.gene.score", "candidate.gene.locus.r2.net", "candidate.gene.locus.r2.net.ensembl", "drug.id.gwas", "drug.num", "drug.id.trial", "drug.gwas.trial.matching", "drug.id.net",
                                                                                  "drug.gwas.net.trial.matching")]
## extend to have each drug as separate entry
drug.table.net <- lapply(1:nrow(drug.table.net), function(x){
  
  ## get entry of interest
  tmp                  <- drug.table.net[x, , drop=F]
  ## get all possible drug IDs (from GWAS and network expansion)
  jj                   <- unique(c(strsplit(tmp$drug.id.gwas, "\\|")[[1]], strsplit(tmp$drug.id.net, "\\|")[[1]]))
  ## report back all possible combinations
  tmp                  <- data.table(tmp, drug.id=jj)
  ## add whether the drug is matching the indication
  jj                   <- unique(c(strsplit(tmp$drug.gwas.trial.matching, "\\|")[[1]], strsplit(tmp$drug.gwas.net.trial.matching, "\\|")[[1]]))
  tmp$match.indication <- tmp$drug.id %in% jj
  ## add number of gwas loci matching
  tmp$num.gwas.loci    <- nrow(res.all[ id == tmp$id[1] & candidate.gene.locus.r2 == tmp$candidate.gene.locus.r2[1]])
  return(tmp)
  
})
## combine
drug.table.net             <- do.call(rbind, drug.table.net)
## n = 73838

## add more information on drugs from OT
drug.table.net             <- merge(drug.table.net, ot.drugs[, c("id", "blackBoxWarning", "description", "drugType", "hasBeenWithdrawn", "isApproved", "name", "synonyms", "tradeNames", "maximumClinicalTrialPhase", 
                                                                 "linkedDiseases.collated", "linkedTargets.collated")],
                                    by.x = "drug.id", by.y = "id")

## indicate the actual overlapping target gene
drug.table.net$target.gene <- apply(drug.table.net[, c("candidate.gene.locus.r2.ensembl", "candidate.gene.locus.r2.net.ensembl", "linkedTargets.collated"), with=F], 1, function(x){
  
  ## report the overlapping gene(s)
  can.gene <- strsplit(x[1], "\\|")[[1]]
  net.gene <- strsplit(x[2], "\\|")[[1]]
  tar.gene <- strsplit(x[3], "\\|")[[1]]
  ## return intersection
  return(paste(intersect(c(can.gene, net.gene), tar.gene), collapse = "|"))
  
})

## indicate whether the target gene is a GWAS locus or a network extension
drug.table.net$type.gene   <- apply(drug.table.net[, c("candidate.gene.locus.r2.ensembl", "target.gene"), with=F], 1, function(x){
  
  ## check for overlap
  jj <- grep(x[1], x[2])
  ## return
  return(length(jj))
  
})

## do in parallel
require(doMC)
registerDoMC(12)

## add information on candidate - network gene relationship
tmp <- mclapply(1:nrow(drug.table.net), function(x){
  
  ## test whether network gene is the one needed
  if(drug.table.net$type.gene[x] == 0){
    
    ## get all available candidate genes at the locus
    gene <- strsplit(drug.table.net$candidate.gene.locus.r2.ensembl[x], "\\|")[[1]]
    ## get all possible target genes
    targ <- strsplit(drug.table.net$target.gene[x], "\\|")[[1]]
    ## get all entries from sig.net referrring to this pair of genes
    tmp <- sig.net[ (ensembl_gene_id.A %in% gene & ensembl_gene_id.B %in% targ) | (ensembl_gene_id.B %in% gene & ensembl_gene_id.A %in% targ), c("EFFECT", "MECHANISM"), with=F]
    ## return
    return(data.table(drug.table.net[x,], net.effect=paste(sort(unique(tmp$EFFECT)), collapse = "|"),
                      net.mechanism=paste(sort(unique(tmp$MECHANISM)), collapse = "|")))
    
    
  }else{
    return(drug.table.net[x,])
  }

}, mc.cores=12)
## combine into one data.table
tmp               <- as.data.table(do.call(plyr::rbind.fill, tmp))

## add symbol for the target gene
gene              <- unique(unlist(lapply(tmp$target.gene, function(x) strsplit(x, "\\|")[[1]])))
## obtain mapping
gene              <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id','external_gene_name'),
                           filters = c('ensembl_gene_id'),
                           values = list(gene),
                           mart = gene.ensembl)
## create new column
tmp$target.symbol <- sapply(tmp$target.gene, function(x){
  ## get all relevant genes
  return(paste(gene$external_gene_name[grep(x, gene$ensembl_gene_id)], collapse = "|"))
})

## re-assign
drug.table.net    <- tmp

################################################
####    examples with matching indications  ####
################################################

#--------------------------------------#
##--  pull indication data from OT  --##
#-------------------------------- -----#

## get all files in the directory
ii <- dir("../downloads/indication/")
ii <- grep("part", ii, value=T)

## run across multiple cores
registerDoMC(10)

## run in parallel
drug.indication <- mclapply(ii, function(x){
  
  ## import each json data set and flatten
  tmp <- ndjson::stream_in(paste0("../downloads/indication/", x))
  
}, mc.cores = 10)
## combine all into one large data frame
drug.indication <- do.call(plyr::rbind.fill, drug.indication)
## keep only what is really needed
drug.indication <- drug.indication[, c("id", "indicationCount", grep("\\.disease|\\.maxPhaseForIndication|\\.efoName", names(drug.indication), value=T))]
## convert to data table
drug.indication <- as.data.table(drug.indication)
## convert from wide to long
drug.indication <- melt(drug.indication, id.vars = c("id", "indicationCount"), measure.vars = list(grep("\\.disease", names(drug.indication), value=T),
                                                                                       grep("\\.maxPhaseForIndication", names(drug.indication), value=T),
                                                                                       grep("\\.efoName", names(drug.indication), value=T)),
                   value.name = c("disease", "maxPhaseForIndication", "efoName"), na.rm = T)
## drop one variable
drug.indication[, variable := NULL]

## add ot id to phecodes
drug.table.net  <- merge(drug.table.net, res.ldsc[, c("id", "ot.id")])

## add this information to drug target table
tmp             <- mclapply(1:nrow(drug.table.net), function(x){
  
  ## get the drug
  d.id <- drug.table.net$drug.id[x]
  ## get all matching disease IDs for the phecode
  p.id <- strsplit(drug.table.net$ot.id[x], "\\|")[[1]]
  ## get all matching entries from indications
  tmp  <- drug.indication[ id == d.id & disease %in% p.id]
  
  ## report findings back
  return(data.table(drug.table.net[x,], maxPhaseForIndication=max(tmp$maxPhaseForIndication), names.indication=paste(tmp$efoName, collapse = "|")))
  
  
}, mc.cores=10)
## combine again
tmp            <- do.call(rbind, tmp)
## replace infinite entries
tmp$maxPhaseForIndication[!is.finite(tmp$maxPhaseForIndication)] <- NA
## reassign
drug.table.net <- tmp

#--------------------------------------#
##--    obtain effect directions    --##
#--------------------------------------#

## identify all drugs with matching indications (market approved)
nrow(drug.table.net[ match.indication == T & !is.na(locus_id.r2) & maxPhaseForIndication == 4])

## create list of variants to be queried to obtain effect sizes
## for indication tested (use maximum phase to sort); drop MHC region!
tmp        <- unique(drug.table.net[ match.indication == T & !is.na(locus_id.r2) & maxPhaseForIndication == 4, c("id", "locus_id.r2")])
## add results from fine-mapping
tmp        <- merge(tmp, res.all[, c("id", "locus_id.r2", "MarkerName", "ID", "CHROM", "GENPOS", "ALLELE0", "ALLELE1", "BETA", "SE", "LOG10P")])
## select only one proxy variant for each locus_id.r2
tmp        <- tmp[order(locus_id.r2, -LOG10P)]
tmp[, ind := 1:.N, by=c("locus_id.r2")]
foo        <- tmp[ ind == 1, c("locus_id.r2", "MarkerName"), with=F]
## add to tmp again
tmp        <- merge(tmp, foo, by=c("locus_id.r2"), suffixes = c(".gwas", ".merged"))

## create sublist of matching indications and associated findings (use common GWAS variants as anchor)
## N.B.: this might drop matchings based purely on gene assignments, but different genetic variants!
drug.matching.effect <- merge(drug.table.net, unique(tmp[, c("locus_id.r2", "MarkerName.merged")]), by="locus_id.r2")
nrow(drug.matching.effect[ match.indication == T & maxPhaseForIndication == 4 ]) ## n = 610

## obtain corresponding effect sizes from variant look-up
var.lookup           <- fread("<path to file>")
## names
names(var.lookup)    <- c("misc", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA")
## create chromosome and phenotype
var.lookup$CHROM     <- as.numeric(gsub(".*:", "", var.lookup$misc)) 
var.lookup$id.phe    <- gsub("/sc-projects/sc-proj-computational-medicine/people/Maik/14_phecode_GWAS/01_GWAS/gwas_results/|\\.allchr.results.gz", "", gsub(":.*", "", var.lookup$misc))
## create MarkerName
var.lookup[, MarkerName := paste0("chr", as.numeric(gsub("X", 23, CHROM)), ":", as.numeric(GENPOS), "_", pmin(ALLELE0,ALLELE1), "_", pmax(ALLELE0,ALLELE1))]

## add to drug matching table
drug.matching.effect <- merge(drug.matching.effect, var.lookup, by.x=c("id", "MarkerName.merged"), by.y=c("id.phe", "MarkerName"))

## delete what is no longer needed
rm(var.lookup); gc()

## create pruned list for repurposing
drug.repurposing     <- drug.matching.effect
## create new indicator of drug -- locus
drug.repurposing[, drug.locus.id := paste0(name, ".", locus_id.r2), ]
## keep only matching example (this may also drop some interesting side effects...)
jj                   <- unique(drug.repurposing[ match.indication == T & maxPhaseForIndication == 4 , "drug.locus.id"]$drug.locus.id)
## subset table accordingly
drug.repurposing     <- drug.repurposing[ drug.locus.id %in% jj]

## drop combinations occurring only once
jj               <- table(drug.repurposing$drug.locus.id)
drug.repurposing <- drug.repurposing[ drug.locus.id %in% names(jj[jj > 1]) ]

## drop combinations with no non-matching
jj               <- do.call(data.frame, aggregate(match.indication ~ drug.locus.id, drug.repurposing, function(x) sum(!x)))
drug.repurposing <- drug.repurposing[ drug.locus.id %in% subset(jj, match.indication > 0)$drug.locus.id ]

#----------------------------------#
##-- create repurposing network --##
#----------------------------------#

## --> edge list 1: locus/gene - disease; include direction <-- ##
edge.rep.gene.disease      <- unique(drug.repurposing[, c("candidate.gene.locus.r2", "phenotype", "BETA")])
## unfold
edge.rep.gene.disease      <- lapply(1:nrow(edge.rep.gene.disease), function(x){
  
  ## get all possible genes by locus
  gene <- strsplit(edge.rep.gene.disease$candidate.gene.locus.r2[x], "\\|")[[1]]
  ## return full table
  return(data.table(candidate.gene.locus.r2=gene, phenotype=edge.rep.gene.disease$phenotype[x], 
                    beta=edge.rep.gene.disease$BETA[x]))
  
})
## combine
edge.rep.gene.disease      <- do.call(rbind, edge.rep.gene.disease)
## n = 165

## generate sign for effect direction
edge.rep.gene.disease$lty  <- ifelse(edge.rep.gene.disease$beta > 0, 1, 2)

## --> edge list 3: drug - gene interaction <-- ##
edge.rep.drug.gene         <- unique(ot.drugs[ linkedTargets.collated != "" & name %in% drug.repurposing$name, c("name", "linkedTargets.collated"), with=F])
## unfold gene names
edge.rep.drug.gene        <- lapply(1:nrow(edge.rep.drug.gene), function(x){
  
  ## unfold all possible genes 
  drug  <- edge.rep.drug.gene$name[x]
  gene  <- strsplit(edge.rep.drug.gene$linkedTargets.collated[x], "\\|")[[1]]
  ## return all combinations
  return(data.frame(drug=drug, gene=gene))
  
})
## combine
edge.rep.drug.gene        <- do.call(rbind, edge.rep.drug.gene) 
## replace Ensembl ID with symbol
tmp.genes                 <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id','external_gene_name', 'gene_biotype'),
                                   filters = c('ensembl_gene_id'),
                                   values = list(unique(edge.rep.drug.gene$gene)),
                                   mart = gene.ensembl)
## add
edge.rep.drug.gene        <- merge(edge.rep.drug.gene, unique(tmp.genes[, c("ensembl_gene_id", "external_gene_name")]), by.x="gene", by.y="ensembl_gene_id")
## reduce to what is really needed
edge.rep.drug.gene        <- edge.rep.drug.gene[, c("drug", "external_gene_name")]
## add dummy type
edge.rep.drug.gene$lty    <- 1
  
## --> edge list 4: drug - disease list (matching indications) <-- ##
edge.rep.drug.phe         <- unique(drug.repurposing[ match.indication == T, c("name", "phenotype")])
## dummy line type
edge.rep.drug.phe$lty     <- 1

## --> edge list 5: gene - gene interaction <-- ##

## reduce to genes included in the network
edge.rep.gene.gene        <- unique(sig.net[ !is.na(external_gene_name.A) & !is.na(external_gene_name.B), c("external_gene_name.A", "external_gene_name.B"), with=F])
## reduce to genes included in the repurposing data
jj                        <- unique(c(edge.rep.locus.gene$candidate.gene.locus.r2, edge.rep.drug.gene$external_gene_name))
edge.rep.gene.gene        <- edge.rep.gene.gene[ external_gene_name.A %in% jj & external_gene_name.B %in% jj]

## define type of the edge
edge.rep.gene.gene$lty    <- 3

## --> combine all edges into one large file <-- ##

## keep only what is really needed and add type of connection
edge.rep.drug.gene$type     <- "drug_target"
edge.rep.drug.phe$type      <- "drug_indication"
edge.rep.gene.gene$type     <- "gene_gene"
edge.rep.gene.disease       <- edge.rep.gene.disease[, c("candidate.gene.locus.r2", "phenotype", "lty")]
edge.rep.gene.disease$type  <- "gene_disease"

## edit names
names(edge.rep.drug.gene) <- names(edge.rep.drug.phe) <- names(edge.rep.gene.gene) <- names(edge.rep.gene.disease) <- c("node1", "node2", "lty", "type")
## combine
edge.repurposing.net        <- rbind(edge.rep.drug.gene, edge.rep.drug.phe, edge.rep.gene.gene, edge.rep.gene.disease)
## drop loops
edge.repurposing.net        <- subset(edge.repurposing.net, node1 != node2)
## n = 832

## create node data 
node.repurposing.net        <- rbind(data.frame(name=unique(c(edge.rep.drug.gene$node1, edge.rep.drug.phe$node1)), type="drug"),
                                     data.frame(name=unique(c(edge.rep.gene.disease$node2, edge.rep.drug.phe$node2)), type="phenotype"),
                                     data.frame(name=unique(c(edge.rep.drug.gene$node2, edge.rep.gene.gene$node1, edge.rep.gene.gene$node2, edge.rep.gene.disease$node1)), type="gene"))
#--------------------------------#
##--   export for webserver   --##
#--------------------------------#

## replace Ensembl ID with symbol
tmp.genes           <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id','external_gene_name', 'gene_biotype'),
                             filters = c('hgnc_symbol'),
                             values = list(unique(subset(node.repurposing.net, type =="gene")$name)),
                             mart = gene.ensembl)

## --> edge data <-- ##
write.table(edge.repurposing.net, "Edges.drug.repurposing.webserver.20230427.txt", sep="\t", row.names=F, quote=F)

## --> drug labels <-- ##
node.drug.webserver <- subset(node.repurposing.net, type == "drug")
## add some information
node.drug.webserver <- merge(node.drug.webserver, unique(ot.drugs[, c("id", "maximumClinicalTrialPhase", "name", "description", "drugType", "blackBoxWarning",
                                                                      "linkedDiseases.collated", "linkedTargets.collated")]), all.x=T, by="name")
## delete one redundant entry
node.drug.webserver <- as.data.table(subset(node.drug.webserver, id != "CHEMBL2105719"))
## add gene names
node.drug.webserver[, linkedTargets.collated.hgnc := sapply(linkedTargets.collated, function(x){
  ## divide
  tmp <- strsplit(x, "\\|")[[1]]
  ## replace (use gene mapping derived below)
  tmp <- sapply(tmp, function(k) tmp.genes$hgnc_symbol[which(tmp.genes$ensembl_gene_id == k)])
  ## return
  return(paste(tmp, collapse = "|"))
})]
## store 
write.table(node.drug.webserver, "Nodes.drug.repurposing.webserver.drugs.20230427.txt", sep="\t", row.names=F, quote=F)

## --> genes <-- ##
node.gene.webserver <- as.data.table(subset(node.repurposing.net, type == "gene"))

## add ensemble ID
node.gene.webserver[, ensembl_gene_id := sapply(node.gene.webserver$name, function(x) paste(unique(tmp.genes$ensembl_gene_id[which(tmp.genes$hgnc_symbol == x)]), collapse = "|"))]
## store
write.table(node.gene.webserver, "Nodes.drug.repurposing.webserver.genes.20230427.txt", sep="\t", row.names=F, quote=F)

## --> phecodes <-- ##
node.phe.webserver  <- as.data.table(subset(node.repurposing.net, type == "phenotype"))
## add some general characteristics
node.phe.webserver  <- merge(node.phe.webserver, res.ldsc[, c("id", "phecode", "phenotype")], by.x = "name", by.y = "phenotype")
## store
write.table(node.phe.webserver, "Nodes.drug.repurposing.webserver.phecode.20230427.txt", sep="\t", row.names=F, quote=F)

################################################
####  network representation of drug table  ####
################################################

#--------------------------------------------#
##--         add effect directions        --##
#--------------------------------------------#

## create list of variants to be queried to obtain effect sizes; drop MHC region and low-confidence assignments
tmp                  <- unique(drug.table.net[ !is.na(locus_id.r2) & (candidate.gene.score >= 1.5 | match.indication == T), c("id", "locus_id.r2")]) ## n = 1735
## add results from fine-mapping
tmp                  <- merge(tmp, res.all[, c("id", "locus_id.r2", "MarkerName", "ID", "CHROM", "GENPOS", "ALLELE0", "ALLELE1", "BETA", "SE", "LOG10P")])
## select only one proxy variant for each locus_id.r2
tmp                  <- tmp[order(locus_id.r2, -LOG10P)]
tmp[, ind := 1:.N, by=c("locus_id.r2")]
foo                  <- tmp[ ind == 1, c("locus_id.r2", "MarkerName"), with=F]
## add to tmp again
tmp                  <- merge(tmp, foo, by=c("locus_id.r2"), suffixes = c(".gwas", ".merged"))

## create sublist of matching indications and associated findings (use common GWAS variants as anchor)
## N.B.: this might drop matchings based purely on gene assignments, but different genetic variants!
drug.all.effect      <- merge(drug.table.net, unique(tmp[, c("locus_id.r2", "MarkerName.merged")]), by="locus_id.r2")

## obtain corresponding effect sizes from variant look-up
var.lookup           <- fread("<path to file>")
## names
names(var.lookup)    <- c("misc", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA")
## create chromosome and phenotype
var.lookup$CHROM     <- as.numeric(gsub(".*:", "", var.lookup$misc)) 
var.lookup$id.phe    <- gsub("/sc-projects/sc-proj-computational-medicine/people/Maik/14_phecode_GWAS/01_GWAS/gwas_results/|\\.allchr.results.gz", "", gsub(":.*", "", var.lookup$misc))
## create MarkerName
var.lookup[, MarkerName := paste0("chr", as.numeric(gsub("X", 23, CHROM)), ":", as.numeric(GENPOS), "_", pmin(ALLELE0,ALLELE1), "_", pmax(ALLELE0,ALLELE1))]

## add to drug matching table
drug.all.effect      <- merge(drug.all.effect, var.lookup, by.x=c("id", "MarkerName.merged"), by.y=c("id.phe", "MarkerName"))

## delete what is no longer needed
rm(var.lookup); gc()


#--------------------------------------------#
##-- generate different edge lists needed --##
#--------------------------------------------#

## --> edge list 1: trait to GWAS gene (may need to be unfolded); restrict to high-confidence assignments <-- ##
edge.trait.gene <- unique(drug.all.effect[, c("candidate.gene.locus.r2", "phenotype", "BETA")])
## unfold
edge.trait.gene <- lapply(1:nrow(edge.trait.gene), function(x){
  
  ## get all possible genes by locus
  gene <- strsplit(edge.trait.gene$candidate.gene.locus.r2[x], "\\|")[[1]]
  ## return full table
  return(data.table(candidate.gene.locus.r2=gene, phenotype=edge.trait.gene$phenotype[x], 
                    beta=edge.trait.gene$BETA[x]))
  
})
## combine
edge.trait.gene       <- do.call(rbind, edge.trait.gene)
## define type of the edge
edge.trait.gene$type  <- "GWAS_gene"
## n = 2044

## generate sign for effect direction
edge.trait.gene$lty   <- ifelse(edge.trait.gene$beta > 0, 1, 2)

## --> edge list 2: gene-gene network (based on Signor); restrict to symbols 'external_gene_name.A' <-- ##
edge.gGene.nGene      <- unique(sig.net[ !is.na(external_gene_name.A) & !is.na(external_gene_name.B), c("external_gene_name.A", "external_gene_name.B"), with=F])
## define type of the edge
edge.gGene.nGene$type <- "network_gene"
## line type
edge.gGene.nGene$lty  <- 3

## --> edge list 3: Drug -> gene target <-- ##
edge.drug.gene        <- unique(ot.drugs[ linkedTargets.collated != "", c("name", "linkedTargets.collated"), with=F])
edge.drug.gene        <- lapply(1:nrow(edge.drug.gene), function(x){
  
  ## unfold all possible genes 
  drug  <- edge.drug.gene$name[x]
  gene  <- strsplit(edge.drug.gene$linkedTargets.collated[x], "\\|")[[1]]
  ## return all combinations
  return(data.frame(drug=drug, gene=gene))
  
})
## combine
edge.drug.gene        <- do.call(rbind, edge.drug.gene) 
## define type of the edge
edge.drug.gene$type   <- "drug_gene"

## replace Ensembl ID with symbol
tmp.genes             <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id','external_gene_name', 'gene_biotype'),
                               filters = c('ensembl_gene_id'),
                               values = list(unique(edge.drug.gene$gene)),
                               mart = gene.ensembl)
## add
edge.drug.gene        <- merge(edge.drug.gene, unique(tmp.genes[, c("ensembl_gene_id", "external_gene_name")]), by.x="gene", by.y="ensembl_gene_id")
## drops some...; reduce to what is really needed
edge.drug.gene        <- edge.drug.gene[, c("drug", "external_gene_name", "type")]
## line type
edge.drug.gene$lty    <- 1

## --> edge list 4: Drug -> trait <-- ##
edge.drug.trait       <- res.ldsc[which(res.ldsc$drug.id != ""), c("phenotype", "drug.id")]
## unfold drugs
edge.drug.trait       <- lapply(1:nrow(edge.drug.trait), function(x){
  
  ## trait
  trait <- edge.drug.trait$phenotype[x]
  ## drug id
  drug  <- strsplit(edge.drug.trait$drug.id[x], "\\|")[[1]]
  ## return
  return(data.frame(drug=drug, trait=trait))
  
})
## combine
edge.drug.trait         <- do.call(rbind, edge.drug.trait)
## replace id with name
edge.drug.trait$drug    <- sapply(edge.drug.trait$drug, function(x) ot.drugs$name[which(ot.drugs$id == x)])
## define type of the edge
edge.drug.trait$type    <- "drug_trait"
## line type
edge.drug.trait$lty     <- 1

#--------------------------------------------#
##--    reduce to what is really needed   --##
#--------------------------------------------#

## drop genes not needed from the Gene - Gene network
edge.gGene.nGene         <- subset(edge.gGene.nGene, external_gene_name.A %in% c(edge.drug.gene$external_gene_name, edge.trait.gene$candidate.gene.locus.r2) & 
                                   external_gene_name.B %in% c(edge.drug.gene$external_gene_name, edge.trait.gene$candidate.gene.locus.r2))

#--------------------------------------------#
##--          generate node data          --##
#--------------------------------------------#

## create common node data set
node.drug.gene.trait      <- data.frame(name=unique(c(edge.drug.gene$external_gene_name, edge.gGene.nGene$external_gene_name.A, edge.gGene.nGene$external_gene_name.B, edge.trait.gene$candidate.gene.locus.r2)), type="gene")
node.drug.gene.trait      <- rbind(node.drug.gene.trait, data.frame(name=unique(c(edge.drug.trait$trait, edge.trait.gene$phenotype)), type="phenotype"))
node.drug.gene.trait      <- rbind(node.drug.gene.trait, data.frame(name=unique(c(edge.drug.trait$drug, edge.drug.gene$drug)), type="drug"))
## add information whether the gene is a GWAS gene
node.drug.gene.trait$gwas <- ifelse(node.drug.gene.trait$name %in% edge.trait.gene$candidate.gene.locus.r2, "gwas_gene", node.drug.gene.trait$type) 

## common names to combine edge files
names(edge.drug.gene)     <- c("node1", "node2", "type", "lty") 
names(edge.drug.trait)    <- c("node1", "node2", "type", "lty")
names(edge.gGene.nGene)   <- c("node1", "node2", "type", "lty")
edge.trait.gene           <- edge.trait.gene[, c("candidate.gene.locus.r2", "phenotype", "type", "lty")]
names(edge.trait.gene)    <- c("node1", "node2", "type", "lty")

## combine into one
edge.drug.gene.trait      <- rbind(edge.drug.gene, edge.drug.trait, edge.gGene.nGene, edge.trait.gene)
## drop loops
edge.drug.gene.trait      <- subset(edge.drug.gene.trait, node1 != node2)
## drop empty nodes
edge.drug.gene.trait      <- subset(edge.drug.gene.trait, node1 != "" & node2 != "")
## drop possible duplications
edge.drug.gene.trait      <- unique(edge.drug.gene.trait)

## export for Cytoscape visualization
write.table(node.drug.gene.trait, "<path to file>", sep="\t", row.names = F, quote = F)
write.table(edge.drug.gene.trait, "<path to file>", sep="\t", row.names = F, quote = F)

##########################################################
####    incorporate QTL data to screen more broadly   ####
##########################################################

#----------------------------#
##-- import coloc results --##
#----------------------------#

## cis-eQTL
res.eqtl <- fread("<path to file>")
## subset to high confidence findings
res.eqtl <- res.eqtl[ PP.H4.abf > .8 & ld.check.top > .8]

## cis-pQTL
res.pqtl <- fread("<path to file>")
## subset to high confidence findings
res.pqtl <- res.pqtl[ PP.H4.abf > .8 & ld.top.overlap > .8]

#----------------------------#
##--   append MR results  --##
#----------------------------#

## package to add Wald ratio estimate
require(MendelianRandomization)
registerDoMC(10)

## --> cis-eQTLs <-- ##
res.eqtl <- mclapply(1:nrow(res.eqtl), function(x){
  
  ## get the relevant data
  mri <- mr_input(bx = res.eqtl$slope[x], bxse = res.eqtl$slope_se[x],
                  by = res.eqtl$Effect.aligned[x], byse = res.eqtl$SE[x])
  ## run the IVW
  mri <- mr_ivw(mri)
  ## return relevant information
  return(data.table(res.eqtl[x,], beta.mr=mri@Estimate, se.mr=mri@StdError, pval.mr=mri@Pvalue))
})
## combine again
res.eqtl <- do.call(rbind, res.eqtl)

## --> cis-pQTLs <-- ##
res.pqtl <- mclapply(1:nrow(res.pqtl), function(x){
  
  ## get the relevant data
  mri <- mr_input(bx = res.pqtl$Effect.protein[x], bxse = res.pqtl$StdErr.protein[x],
                  by = res.pqtl$Effect.trait[x], byse = res.pqtl$StdErr.trait[x])
  ## run the IVW
  mri <- mr_ivw(mri)
  ## return relevant information
  return(data.table(res.pqtl[x,], beta.mr=mri@Estimate, se.mr=mri@StdError, pval.mr=mri@Pvalue))
})
## combine again
res.pqtl <- do.call(rbind, res.pqtl)

#---------------------------------#
##--  add phecode information  --##
#---------------------------------#

## add druggable information mapped to phecodes
tmp.ldsc <- fread("<path to file>")
res.ldsc <- merge(res.ldsc, tmp.ldsc, by=intersect(names(res.ldsc), names(tmp.ldsc)))
  
## cis-eQTL
res.eqtl <- merge(res.eqtl, res.ldsc[, c("id", "category", "phecode", "phenotype", "cl", "p.care", "srt.category", "ot.id", "drug.id", "drug.phase")])
## cis-pQTL
res.pqtl <- merge(res.pqtl, res.ldsc[, c("id", "category", "phecode", "phenotype", "cl", "p.care", "srt.category", "ot.id", "drug.id", "drug.phase")],
                  by=intersect(names(res.pqtl), c("id", "category", "phecode", "phenotype", "cl", "p.care", "srt.category", "ot.id", "drug.id", "drug.phase")))

#---------------------------------#
##--    map drug information   --##
#---------------------------------#

## create new target list
tmp.drugs <- ot.drugs[ linkedTargets.collated != ""]
## drop targets that have already been withdrawn
tmp.drugs <- tmp.drugs[ hasBeenWithdrawn == F]

## unfold by unqiue targets
tmp.drugs <- mclapply(1:nrow(tmp.drugs), function(x){
  
  ## get all relevant entries
  jj <- strsplit(tmp.drugs$linkedTargets.collated[x], "\\|")[[1]]
  ## return
  return(data.table(tmp.drugs[x,], ensembl_gene_id = jj))
}, mc.cores = 10)
## combine again
tmp.drugs <- do.call(rbind, tmp.drugs)

## --> cis-eQTL <-- ##
res.eqtl$drug.id.target <- sapply(res.eqtl$ensembl_gene_id, function(x){
  ## get all relevant durgs (if any)
  tmp <- tmp.drugs[ ensembl_gene_id == x]
  return(paste(tmp$id, collapse = "|"))
})
## subset to those with at least one mapping
res.eqtl                <- res.eqtl[ drug.id.target != ""]

## --> cis-pQTL <-- ##
res.pqtl$drug.id.target <- sapply(res.pqtl$ensembl_gene_id, function(x){
  ## get all relevant durgs (if any)
  tmp <- tmp.drugs[ ensembl_gene_id == x]
  return(paste(tmp$id, collapse = "|"))
})
## subset to those with at least one mapping
res.pqtl                <- res.pqtl[ drug.id.target != ""]

#---------------------------------#
##--    add locus information  --##
#---------------------------------#

## add locus_id.r2 information to coloc results to identify pleiotropic loci
res.eqtl  <- merge(res.eqtl, res.all[, c("id", "MarkerName", "locus_id.r2")], by.x = c("id", "MarkerName.phe"), by.y = c("id", "MarkerName"))
res.pqtl  <- merge(res.pqtl, res.all[, c("id", "MarkerName", "locus_id.r2")], by.x = c("id", "MarkerName.phe"), by.y = c("id", "MarkerName"))

## --> reduce cis-eQTLs <-- ##

## keep only one for each id - gene - tissue pair
res.eqtl  <- res.eqtl[ order(id, ensembl_gene_id, tissue, -PP.H4.abf)]
res.eqtl[, ind := 1:.N, by=c("id", "ensembl_gene_id", "tissue")]
res.eqtl  <- res.eqtl[ ind == 1]

## drop loci with more than one gene
loci.eqtl <- do.call(data.table, aggregate(ensembl_gene_id ~ locus_id.r2, res.eqtl, function(x) paste(sort(unique(x)), collapse = "|")))
## compress again by shared genes
loci.eqtl <- do.call(data.table, aggregate(locus_id.r2 ~ ensembl_gene_id, loci.eqtl, function(x) paste(sort(unique(x)), collapse = "|")))

## --> reduce cis-eQTLs <-- ##

## keep only one for each id - gene - tissue pair
res.pqtl  <- res.pqtl[ order(id, pheno, -PP.H4.abf)]
res.pqtl[, ind := 1:.N, by=c("id", "pheno")]
res.pqtl  <- res.pqtl[ ind == 1]

## drop loci with more than one gene
loci.pqtl <- do.call(data.table, aggregate(ensembl_gene_id ~ locus_id.r2, res.pqtl, function(x) paste(sort(unique(x)), collapse = "|")))
## compress again by shared genes
loci.pqtl <- do.call(data.table, aggregate(locus_id.r2 ~ ensembl_gene_id, loci.pqtl, function(x) paste(sort(unique(x)), collapse = "|")))

#---------------------------------#
##--  add matching indication  --##
#---------------------------------#

## --> cis-eQTLs <-- ##
res.eqtl[, match.indication := apply(res.eqtl[, c("drug.id", "drug.id.target"), with=F], 1, function(x){
  ## get all drugs
  ind.d <- strsplit(x[1], "\\|")[[1]]
  gen.d <- strsplit(x[2], "\\|")[[1]]
  ## report back
  return(ifelse(length(intersect(ind.d, gen.d)) > 0, T, F))
})]

## --> cis-pQTLs <-- ##
res.pqtl[, match.indication := apply(res.pqtl[, c("drug.id", "drug.id.target"), with=F], 1, function(x){
  ## get all drugs
  ind.d <- strsplit(x[1], "\\|")[[1]]
  gen.d <- strsplit(x[2], "\\|")[[1]]
  ## report back
  return(ifelse(length(intersect(ind.d, gen.d)) > 0, T, F))
})]

#-------------------------------------#
##-- unfold and add mode of action --##
#-------------------------------------#

## unfold by mapping drugs - cis-eQTLs
res.eqtl <- lapply(1:nrow(res.eqtl), function(x){
  
  ## get all matching drugs
  drug <- strsplit(res.eqtl$drug.id.target[x], "\\|")[[1]]
  ## return
  return(data.table(res.eqtl[x,], id.drug = drug))
  
})
## combine
res.eqtl <- do.call(rbind, res.eqtl)

## unfold by mapping drugs - cis-pQTLs
res.pqtl <- lapply(1:nrow(res.pqtl), function(x){
  
  ## get all matching drugs
  drug <- strsplit(res.pqtl$drug.id.target[x], "\\|")[[1]]
  ## return
  return(data.table(res.pqtl[x,], id.drug = drug))
  
})
## combine
res.pqtl <- do.call(rbind, res.pqtl)

## --> add mechanism of action <-- ##

## cis-eQTLs
res.eqtl <- merge(res.eqtl, ot.drugs.moa, by.x="id.drug", by.y="id", all.x=T)
res.eqtl <- merge(res.eqtl, ot.drugs[, c("id", "blackBoxWarning", "description", "name", "linkedDiseases.collated", "linkedTargets.collated")], by.x = "id.drug", by.y = "id")
## careful, not all have a mode of action, and some have multiple

## cis-pQTLs
res.pqtl <- merge(res.pqtl, ot.drugs.moa, by.x="id.drug", by.y="id", all.x=T)
res.pqtl <- merge(res.pqtl, ot.drugs[, c("id", "blackBoxWarning", "description", "name", "linkedDiseases.collated", "linkedTargets.collated")], by.x = "id.drug", by.y = "id")
## careful, not all have a mode of action, and some have multiple

## --> ensure the MoA fits to the right target <-- ##

## cis-eQTLs
res.eqtl[, fit.moa := apply(res.eqtl[, c("ensembl_gene_id", "targets"), with=F], 1, function(x){
  ## get target genes for the MoA
  m.gen <- strsplit(x[2], "\\|")[[1]]
  ## return whether gene is included
  return(x[1] %in% m.gen)
})]
## drop none matching ones
res.eqtl <- res.eqtl[ fit.moa == T]

## cis-eQTLs
res.pqtl[, fit.moa := apply(res.pqtl[, c("ensembl_gene_id", "targets"), with=F], 1, function(x){
  ## get target genes for the MoA
  m.gen <- strsplit(x[2], "\\|")[[1]]
  ## return whether gene is included
  return(x[1] %in% m.gen)
})]
## drop none matching ones
res.pqtl <- res.pqtl[ fit.moa == T]

#-------------------------------------#
##--  compare MR finding with MOA  --##
#-------------------------------------#

## create data frame to do so
tmp.moa  <- data.table(actionType = c("ACTIVATOR", "AGONIST", "ANTAGONIST", "ANTISENSE INHIBITOR", "BINDING AGENT", "BLOCKER", 
                                      "HYDROLYTIC ENZYME", "INHIBITOR", "INVERSE AGONIST", "MODULATOR", "NEGATIVE ALLOSTERIC MODULATOR",
                                      "OPENER", "OTHER", "PARTIAL AGONIST", "POSITIVE ALLOSTERIC MODULATOR", "POSITIVE MODULATOR", "PROTEOLYTIC ENZYME",
                                      "RNAI INHIBITOR", "STABILISER", "VACCINE ANTIGEN"),
                       effect.on.target=c(1, 1, -1, -1,  -1, -1, -1, -1, -1, 1, -1, 1, 0, 1, 1, 1, -1, -1, 1, -1))

## add to cis-eQTL
res.eqtl <- merge(res.eqtl, tmp.moa, by="actionType")
## add to cis-pQTL
res.pqtl <- merge(res.pqtl, tmp.moa, by="actionType")
  
## --> compare association directions and drug effects <-- ##

## cis-eQTL
res.eqtl[, type.effect := ifelse(effect.on.target == 0, "unkown", ifelse(effect.on.target == sign(beta.mr), "adverse", "beneficial"))]
## cis-pQTL
res.pqtl[, type.effect := ifelse(effect.on.target == 0, "unkown", ifelse(effect.on.target == sign(beta.mr), "adverse", "beneficial"))]

#-------------------------------------#
##-- recompute indication matching --##
#-------------------------------------#

## cis-eQTLs
res.eqtl[, match.indication := apply(res.eqtl[, c("id.drug", "drug.id"), with=F], 1, function(x){
  ## get all drugs
  ind.d <- strsplit(x[2], "\\|")[[1]]
  ## report back
  return(x[1] %in% ind.d)
})]

## cis-pQTLs
res.pqtl[, match.indication := apply(res.pqtl[, c("id.drug", "drug.id"), with=F], 1, function(x){
  ## get all drugs
  ind.d <- strsplit(x[2], "\\|")[[1]]
  ## report back
  return(x[1] %in% ind.d)
})]

## add whether or not the disease is mostly seen in primary care
res.eqtl <- merge(res.eqtl, res.ldsc[, c("id", "case.ratio.pcare")], by="id")
res.pqtl <- merge(res.pqtl, res.ldsc[, c("id", "case.ratio.pcare")], by="id")

#-------------------------------------#
##--        compress by drug       --##
#-------------------------------------#

## compress by drug - gene combination; cis-eQTL
drug.gene.eqtl <- unique(res.eqtl[, c("id.drug", "ensembl_gene_id", "gene_name", "name", "actionType", "mechanismOfAction", "targetName", "targetType", "targets", "blackBoxWarning",
                                     "description", "effect.on.target", "linkedDiseases.collated")])

## add relevant columns
drug.gene.eqtl <- lapply(1:nrow(drug.gene.eqtl), function(x){
  
  ## get the ids to grep all relevant entries
  drug  <- drug.gene.eqtl$id.drug[x]
  gene  <- drug.gene.eqtl$ensembl_gene_id[x]
  tar   <- drug.gene.eqtl$targets[x]
  ## get all relevant entries
  tmp   <- res.eqtl[ id.drug == drug & ensembl_gene_id == gene & targets == tar]
  ## devide into beneficial and adverse effects
  tmp.b <- tmp[ type.effect == "beneficial"]
  tmp.a <- tmp[ type.effect == "adverse"]
  tmp.m <- tmp[match.indication == T]
  ## generate output file
  return(data.table(drug.gene.eqtl[x,], 
                    beneficial.phecodes.tissue=paste0(apply(tmp.b[, c("phenotype", "id", "tissue"), with=F], 1, function(k) paste(k, collapse = "-")), collapse = "|"), 
                    beneficial.phecodes=paste(sort(unique(tmp.b$phenotype)), collapse = "|"),
                    adverse.phecodes.tissue=paste0(apply(tmp.a[, c("phenotype", "id", "tissue"), with=F], 1, function(k) paste(k, collapse = "-")), collapse = "|"),
                    adverse.phecodes=paste(sort(unique(tmp.a$phenotype)), collapse = "|"),
                    match.indication=paste(sort(unique(tmp.m$phenotype)), collapse = "|")))
})
## combine again
drug.gene.eqtl <- do.call(rbind, drug.gene.eqtl)

## --> now for cis-pQTLs <-- ##
drug.gene.pqtl <- unique(res.pqtl[, c("id.drug", "ensembl_gene_id", "hgnc_symbol", "name", "actionType", "mechanismOfAction", "targetName", "targetType", "targets", "blackBoxWarning",
                                     "description", "effect.on.target", "linkedDiseases.collated")])

## add relevant columns
drug.gene.pqtl <- lapply(1:nrow(drug.gene.pqtl), function(x){
  
  ## get the ids to grep all relevant entries
  drug  <- drug.gene.pqtl$id.drug[x]
  gene  <- drug.gene.pqtl$ensembl_gene_id[x]
  tar   <- drug.gene.pqtl$targets[x]
  ## get all relevant entries
  tmp   <- res.pqtl[ id.drug == drug & ensembl_gene_id == gene & targets == tar]
  ## devide into beneficial and adverse effects
  tmp.b <- tmp[ type.effect == "beneficial"]
  tmp.a <- tmp[ type.effect == "adverse"]
  tmp.m <- tmp[match.indication == T]
  ## generate output file
  return(data.table(drug.gene.pqtl[x,], 
                    beneficial.phecodes.tissue=paste0(apply(tmp.b[, c("phenotype", "id"), with=F], 1, function(k) paste(k, collapse = "-")), collapse = "|"), 
                    beneficial.phecodes=paste(sort(unique(tmp.b$phenotype)), collapse = "|"),
                    adverse.phecodes.tissue=paste0(apply(tmp.a[, c("phenotype", "id"), with=F], 1, function(k) paste(k, collapse = "-")), collapse = "|"),
                    adverse.phecodes=paste(sort(unique(tmp.a$phenotype)), collapse = "|"),
                    match.indication=paste(sort(unique(tmp.m$phenotype)), collapse = "|")))
})
## combine again
drug.gene.pqtl <- do.call(rbind, drug.gene.pqtl)

## --> combine into one large file <-- ##

## combine both resources
names(drug.gene.pqtl) <- gsub("hgnc_symbol", "gene_name", names(drug.gene.pqtl))
drug.gene.QTL         <- merge(drug.gene.eqtl, drug.gene.pqtl, by=c("id.drug", "ensembl_gene_id", "gene_name", "name", "actionType", "mechanismOfAction", "targetName", "targetType", "targets", "blackBoxWarning",
                                                                    "description", "effect.on.target", "linkedDiseases.collated"),
                               suffixes = c(".eqtl", ".pqtl"), all=T, )
## fuse beneficial effects
drug.gene.QTL[, beneficial.phecodes := apply(drug.gene.QTL[, c("beneficial.phecodes.eqtl", "beneficial.phecodes.pqtl")], 1, function(x){
  ## split up
  if(!is.na(x[1])){
    ii <- strsplit(x[1], "\\|")[[1]]
  }else{
    ii <- ""
  }
  if(!is.na(x[2])){
    jj <- strsplit(x[2], "\\|")[[1]]
  }else{
    jj <- ""
  }
  ## return
  return(paste(sort(unique(c(ii,jj))), collapse = "|"))
})]

## fuse adverse effects
drug.gene.QTL[, adverse.phecodes := apply(drug.gene.QTL[, c("adverse.phecodes.eqtl", "adverse.phecodes.pqtl")], 1, function(x){
  ## split up
  if(!is.na(x[1])){
    ii <- strsplit(x[1], "\\|")[[1]]
  }else{
    ii <- ""
  }
  if(!is.na(x[2])){
    jj <- strsplit(x[2], "\\|")[[1]]
  }else{
    jj <- ""
  }
  ## return
  return(paste(sort(unique(c(ii,jj))), collapse = "|"))
})]

## tissue dependent effects
drug.gene.QTL[, tissue.diff := apply(drug.gene.QTL[, c("beneficial.phecodes", "adverse.phecodes"), with=F], 1, function(x){
  ## get all beneficial effects
  ii <- strsplit(x[1], "\\|")[[1]]
  ## get all adverse effects
  jj <- strsplit(x[2], "\\|")[[1]]
  ## test whether there is any overlap
  jj <- intersect(ii, jj)
  return(length(jj))
})]

## --> write to file to create supplemental tables <-- ##

## based from a drugs perspective
write.table(drug.gene.QTL, "Putative.drug.repurposing.adverse.effects.QTL.fused.20230103.txt", sep="\t", row.names=F)

## collate QTL evidence
tmp.eqtl        <- res.eqtl[, c("id", "category", "phenotype", "name", "type.effect", "ensembl_gene_id", "gene_name", "tissue", "match.indication", "mechanismOfAction", 
                                "blackBoxWarning", "description", "actionType", "id.drug", "MarkerName.eQTL", "ID", "locus_id.r2", "Allele1", "Allele2", "Effect.aligned", "SE", "LOG10P",
                                "slope", "slope_se", "pval_nominal", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf", "beta.mr", "se.mr", "pval.mr")]
tmp.pqtl        <- res.pqtl[, c("id", "category", "phenotype", "name", "type.effect", "ensembl_gene_id", "hgnc_symbol", "match.indication", "mechanismOfAction", 
                                "blackBoxWarning", "description", "actionType", "id.drug", "MarkerName.pQTL", "rsid", "locus_id.r2", "ALLELE0", "ALLELE1", "Effect.trait", "StdErr.trait", "log10p.trait",
                                "Effect.protein", "StdErr.protein", "Pvalue.protein", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf", "beta.mr", "se.mr", "pval.mr")]
## edit names
names(tmp.eqtl) <- c("id", "category", "phenotype", "name", "type.effect", "ensembl_gene_id", "gene_name", "tissue", "match.indication", "mechanismOfAction", 
                     "blackBoxWarning", "description", "actionType", "id.drug", "MarkerName", "rsid", "locus_id.r2", "non_effect_allele", "effect_allele", "Effect.trait", "SE.trait", "log10p.trait",
                     "Effect.qtl", "StdErr.qtl", "Pvalue.qtl", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf", "beta.mr", "se.mr", "pval.mr")
names(tmp.pqtl) <- c("id", "category", "phenotype", "name", "type.effect", "ensembl_gene_id", "gene_name", "match.indication", "mechanismOfAction", 
                     "blackBoxWarning", "description", "actionType", "id.drug", "MarkerName", "rsid", "locus_id.r2", "non_effect_allele", "effect_allele", "Effect.trait", "SE.trait", "log10p.trait",
                     "Effect.qtl", "StdErr.qtl", "Pvalue.qtl", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf", "beta.mr", "se.mr", "pval.mr")
## add some columns
tmp.eqtl[, type.qtl := "cis_eQTL"]
tmp.pqtl[, tissue := "Blood"]
tmp.pqtl[, type.qtl := "cis_pQTL"]

## combine
res.qtl         <- rbind(tmp.eqtl, tmp.pqtl)

## add OR column for MR results
res.qtl[, OR.mr := paste0(sprintf("%.2f", exp(beta.mr)), " (", sprintf("%.2f", exp(beta.mr - 1.96 * se.mr)), ";", sprintf("%.2f", exp(beta.mr + 1.96 * se.mr)),")")]

## write to file
write.table(res.qtl, "<path to file>", row.names = F, sep="\t")

#####################################################################
#####################################################################
####                  Revision NatGen 05/10/2023                 ####
#####################################################################
#####################################################################

#########################################
####   enrichment of druggable genes ####
#########################################

## list of druggable genes
drug.genes <- unique(unlist(lapply(ot.drugs$linkedTargets.collated, function(x) strsplit(x, "\\|")[[1]]))) ## N = 1462

## test whether druggable genes are enriched in different 
## bins of gene assignment
drug.enrich <- lapply(1:4, function(x){
  
  ## get all genes assigned within this bin
  jj <- res.all[ candidate.gene.score > x-1 & candidate.gene.score <= ifelse(x == 4, 8, x) ]$candidate.gene.locus.r2.ensembl
  jj <- unique(unlist(lapply(jj, function(x) strsplit(x, "\\|")[[1]])))
  
  ## compute enrichment for druggable genes
  
  ## druggable and in phecode data
  d1  <- length(intersect(jj, drug.genes))
  ## druggable and not in phecode data
  d2  <- length(drug.genes[!(drug.genes %in% jj)])
  ## not druggable and in phecode data
  d3  <- length(jj[!(jj %in% drug.genes)])
  ## not druggable and not in phecode data (19,817  Genome Res. 22, 1760–1774 (2012))
  d4  <- 19817 - length(unique(c(jj, drug.genes)))
  
  ## test for enrichment
  print(matrix(c(d1, d2, d3, d4), 2, 2, byrow = T))
  enr <- fisher.test(matrix(c(d1, d2, d3, d4), 2, 2, byrow = T))
  
  ## return information needed
  return(data.frame(score.bin=x, or=enr$estimate, cil=enr$conf.int[1], ciu=enr$conf.int[2], pval=enr$p.value))
  
  
})
## combine everything
drug.enrich <- do.call(rbind, drug.enrich)

## --> simple figure <-- ##

#########################################
####     fraction realted diseases   ####
#########################################

## import disease network
disease.net      <- fread("/sc-projects/sc-proj-computational-medicine/people/Maik/14_phecode_GWAS/15_COVID_prediction/input/UKB.phecode.partial.correlation.network.20230926.txt")

## create column, whether the index conditions is related (same ICD-10 chapter) to the 'new' indication
tmp              <- unique(drug.repurposing[ match.indication == T, c("name", "phenotype")]) ## still multiple indications by drug
## add category
tmp              <- merge(tmp, res.ldsc[, c("phenotype", "category")])
## add network diseases
tmp[, network.disease := sapply(tmp$phenotype, function(x){
  ## get subnetwork
  jj <- disease.net[ phenotype.1 == x | phenotype.2 == x]
  ## return
  return(paste(sort(unique(c(jj$phenotype.1, jj$phenotype.2))), collapse = "|"))
})]

## add disease category to repurposing information
drug.repurposing <- merge(drug.repurposing, res.ldsc[, c("phenotype", "category")], by="phenotype")

## add additional information on index disease
drug.repurposing[, related.to.index.condition := apply(drug.repurposing[, c("name", "phenotype", "category")], 1, function(x){
  
  ## get relevant index conditions
  jj <- subset(tmp, name == x[1])

  ## test whether it is related
  if(x[2] %in% jj$phenotype | x[3] %in% jj$category | x[2] %in% unlist(lapply(jj$network.disease, function(c) strsplit(c, "\\|")[[1]]))){
    return("yes")
  }else{
    return("no")
  }
  
})]

## table
table(drug.repurposing$related.to.index.condition, drug.repurposing$match.indication)
#     FALSE TRUE
# no    138    0
# yes   253  212

##########################################
#### revise list of molecular targets ####
##########################################

## flag loci for which the coloc candidate gene is different
res.qtl     <- merge(res.qtl, unique(res.all[, c("locus_id.r2", "candidate.gene.locus.r2", "candidate.gene.locus.r2.ensembl")]))
## create flag
res.qtl[, locus.gene.map := apply(res.qtl[, c("gene_name", "candidate.gene.locus.r2"), with=F], 1, function(x){
  ## get all possible candidate genes
  jj <- strsplit(x[2], "\\|")[[1]]
  ## return results
  return(ifelse(x[1] %in% jj, "yes", "no"))
})]
## table
table(res.qtl$locus.gene.map)
#  no  yes 
# 961 3096 

## --> look at PubMed <-- ##
set_entrez_key("be6f0336289b95f13fae71ae52b40ae52e08")

## flag whether phenotype and medication have been commonly mentioned in pubmed abtracts
tmp.pubmed  <- unique(res.qtl[, c("name", "phenotype")])
## n = 1572 combinations
tmp.pubmed  <- lapply(1:nrow(tmp.pubmed), function(x){
  
  ## define query; restrict to abstract or title
  qy       <- paste("(", tmp.pubmed$name[x], ") AND (", tmp.pubmed$phenotype[x],")")
  ## print search
  print(qy)
  ## perform the search
  r_search <- entrez_search(db="pubmed", term=qy, retmax = 10)
  print(r_search$QueryTranslation)
  print(r_search$count)
  ## return output
  return(data.table(tmp.pubmed[x,], pubmed.entries=r_search$count))
  
})
## combine
tmp.pubmed  <- do.call(rbind, tmp.pubmed)

## add to molecular targets
res.qtl     <- merge(as.data.table(res.qtl), as.data.table(tmp.pubmed), by=c("name", "phenotype"))

## --> drop drug - target mappings w/o evidence from gene assignment <-- ##

## get list of gene/drug pairings to keep (keep only positive examples or those with consistent evidence of causal gene annotation)
tmp.keep    <- unique(res.qtl[locus.gene.map == "yes" , c("name", "gene_name")])
## n = 715 unique drug - gene combinations

## subset the list accordingly
# res.qtl.rev <- res.qtl[ paste(name, gene_name, sep="_") %in% apply(tmp.keep, 1, paste, collapse="_")]
res.qtl.rev <- res.qtl[ locus.gene.map == "yes" ]
## n = 3096

## how many unique genes
length(unique(res.qtl.rev$gene_name)) ## n = 105
## by molecular entity
aggregate(gene_name ~ type.qtl, res.qtl.rev, function(x) length(unique(x)))
## how many diseases
length(unique(res.qtl.rev$phenotype)) ## n = 99

## how many unique drugs
length(unique(res.qtl.rev$name)) ## n = 685
## how many unique genes

## write to file
write.table(res.qtl.rev, "Results.fused.QTL.coloc.drugs.20231109.txt", row.names = F, sep="\t")

#-------------------------------------#
##--        compress by drug       --##
#-------------------------------------#

## compress by drug - gene combination; cis-eQTL
drug.gene.QTL.rev <- unique(res.qtl.rev[, c("id.drug", "ensembl_gene_id", "gene_name", "name", "actionType", "mechanismOfAction", "blackBoxWarning", "description")])
## n = 720 

## add locus map, to ensure that no unspecific effects are picked up
res.pqtl          <- merge(res.pqtl, unique(res.all[, c("locus_id.r2", "candidate.gene.locus.r2", "candidate.gene.locus.r2.ensembl")]), by = "locus_id.r2")
## create flag
res.pqtl[, locus.gene.map := apply(res.pqtl[, c("hgnc_symbol", "candidate.gene.locus.r2"), with=F], 1, function(x){
  ## get all possible candidate genes
  jj <- strsplit(x[2], "\\|")[[1]]
  ## return results
  return(ifelse(x[1] %in% jj, "yes", "no"))
})]

## add locus map, to ensure that no unspecific effects are picked up
res.eqtl          <- merge(res.eqtl, unique(res.all[, c("locus_id.r2", "candidate.gene.locus.r2", "candidate.gene.locus.r2.ensembl")]), by = "locus_id.r2")
## create flag
res.eqtl[, locus.gene.map := apply(res.eqtl[, c("gene_name", "candidate.gene.locus.r2"), with=F], 1, function(x){
  ## get all possible candidate genes
  jj <- strsplit(x[2], "\\|")[[1]]
  ## return results
  return(ifelse(x[1] %in% jj, "yes", "no"))
})]

## add relevant columns
drug.gene.QTL.rev <- lapply(1:nrow(drug.gene.QTL.rev), function(x){
  
  ## get the ids to grep all relevant entries
  drug  <- drug.gene.QTL.rev$id.drug[x]
  gene  <- drug.gene.QTL.rev$ensembl_gene_id[x]
  
  ## get all relevant cis-eQTL entries
  e.tmp   <- res.eqtl[ id.drug == drug & ensembl_gene_id == gene & locus.gene.map == "yes" ]
  ## devide into beneficial and adverse effects
  e.tmp.b <- e.tmp[ type.effect == "beneficial"]
  e.tmp.a <- e.tmp[ type.effect == "adverse"]
  e.tmp.m <- e.tmp[match.indication == T]
  
  ## get all relevant cis-pQTL entries
  p.tmp   <- res.pqtl[ id.drug == drug & ensembl_gene_id == gene & locus.gene.map == "yes" ]
  ## devide into beneficial and adverse effects
  p.tmp.b <- p.tmp[ type.effect == "beneficial"]
  p.tmp.a <- p.tmp[ type.effect == "adverse"]
  p.tmp.m <- p.tmp[match.indication == T]
  
  ## generate output file
  return(data.table(drug.gene.QTL.rev[x,], 
                    beneficial.phecodes.tissue.eqtl=paste0(apply(e.tmp.b[, c("phenotype", "id", "tissue"), with=F], 1, function(k) paste(k, collapse = "-")), collapse = "|"),
                    beneficial.phecodes.tissue.pqtl=paste0(apply(p.tmp.b[, c("phenotype", "id"), with=F], 1, function(k) paste(k, collapse = "-")), collapse = "|"),
                    beneficial.phecodes.eqtl=paste(sort(unique(e.tmp.b$phenotype)), collapse = "|"),
                    beneficial.phecodes.pqtl=paste(sort(unique(p.tmp.b$phenotype)), collapse = "|"),
                    adverse.phecodes.tissue.eqtl=paste0(apply(e.tmp.a[, c("phenotype", "id", "tissue"), with=F], 1, function(k) paste(k, collapse = "-")), collapse = "|"),
                    adverse.phecodes.tissue.pqtl=paste0(apply(p.tmp.a[, c("phenotype", "id"), with=F], 1, function(k) paste(k, collapse = "-")), collapse = "|"),
                    adverse.phecodes.eqtl=paste(sort(unique(e.tmp.a$phenotype)), collapse = "|"),
                    adverse.phecodes.pqtl=paste(sort(unique(p.tmp.a$phenotype)), collapse = "|"),
                    match.indication.eqtl=paste(sort(unique(e.tmp.m$phenotype)), collapse = "|"),
                    match.indication.pqtl=paste(sort(unique(p.tmp.m$phenotype)), collapse = "|")))
})
## combine again
drug.gene.QTL.rev <- do.call(rbind, drug.gene.QTL.rev)
## n = 720

## how many unique drugs
length(unique(drug.gene.QTL.rev$name)) ## n = 677
## how many unique genes
length(unique(drug.gene.QTL.rev$gene_name)) ## n = 105
## how many unique phecodes
length(unique(res.qtl.rev$phenotype))


## fuse beneficial effects
drug.gene.QTL.rev[, beneficial.phecodes := apply(drug.gene.QTL.rev[, c("beneficial.phecodes.eqtl", "beneficial.phecodes.pqtl")], 1, function(x){
  ## split up
  if(!is.na(x[1])){
    ii <- strsplit(x[1], "\\|")[[1]]
  }else{
    ii <- ""
  }
  if(!is.na(x[2])){
    jj <- strsplit(x[2], "\\|")[[1]]
  }else{
    jj <- ""
  }
  ## return
  return(paste(sort(unique(c(ii,jj))), collapse = "|"))
})]

## fuse adverse effects
drug.gene.QTL.rev[, adverse.phecodes := apply(drug.gene.QTL.rev[, c("adverse.phecodes.eqtl", "adverse.phecodes.pqtl")], 1, function(x){
  ## split up
  if(!is.na(x[1])){
    ii <- strsplit(x[1], "\\|")[[1]]
  }else{
    ii <- ""
  }
  if(!is.na(x[2])){
    jj <- strsplit(x[2], "\\|")[[1]]
  }else{
    jj <- ""
  }
  ## return
  return(paste(sort(unique(c(ii,jj))), collapse = "|"))
})]

## tissue dependent effects
drug.gene.QTL.rev[, tissue.diff := apply(drug.gene.QTL.rev[, c("beneficial.phecodes", "adverse.phecodes"), with=F], 1, function(x){
  ## get all beneficial effects
  ii <- strsplit(x[1], "\\|")[[1]]
  ## get all adverse effects
  jj <- strsplit(x[2], "\\|")[[1]]
  ## test whether there is any overlap
  jj <- intersect(ii, jj)
  return(length(jj))
})]

## create common matching category
drug.gene.QTL.rev[, matching.indication := apply(drug.gene.QTL.rev[, c("match.indication.eqtl", "match.indication.pqtl"), with=F], 1, function(x){
  x <- unique(x)
  x <- x[x != ""]
  return(paste(unique(x), collapse = "|"))
})]

## flag whether there is a misalignment of effect directions for indications
drug.gene.QTL.rev[, effect_conflict := apply(drug.gene.QTL.rev[, c("matching.indication", "adverse.phecodes"), with=F], 1, function(x){
  ## get all beneficial effects
  ii <- strsplit(x[1], "\\|")[[1]]
  ## get all adverse effects
  jj <- strsplit(x[2], "\\|")[[1]]
  ## test whether there is any overlap
  jj <- intersect(ii, jj)
  if(length(jj) > 0){
    return("yes")
  }else{
    return("no")
  }
})]
## have a look
table(drug.gene.QTL.rev$effect_conflict)
#  no yes 
# 674  46 

## --> add reported side effects!! <-- ##

## add reported side effects from Chris Finan (https://cfinan.gitlab.io/bio-misc/scripts/drug_lookups.html#bm-drug-target-effects)
ade.summary   <- fread("drug_target_pairs_all_effects_20231020.txt")
## only 621 mappings by compound name

## query drug effects in bioportal to obtain possible mappings
source("../scripts/search_bioportal.R")
## api key
api_key       <- "a745bbc9-9640-496d-958b-c8929008da90"
## get side effects, [sample(1:length(ade.summary$drug_effect), 20)]
ade.mapping   <- lapply(unique(ade.summary$drug_effect), function(x){
  print(x)
  ## get the results from the web query
  tmp     <- search_bioportal(x, "ICD10,SNOMEDCT,MEDDRA", api_key)
  ## add the ADE of interest
  tmp$ade <- x
  ## return
  return(tmp)
})
## clean errors
ade.mapping   <- Filter(function(x) is.data.frame(x), ade.mapping)
## combine everything
ade.mapping   <- do.call(plyr::rbind.fill, ade.mapping)
## drop non-matching ones
ade.mapping   <- ade.mapping[ !is.na(source)]
## keep only what is needed
ade.mapping   <- unique(ade.mapping[, c("notation", "prefLabel", "source", "ade")])

## import possible mapping tables
meddra.icd10  <- fread("/sc-projects/sc-proj-computational-medicine/people/Maik/14_phecode_GWAS/19_phecode_disease_mapping/input/Mapping.MEDDRA.ICD10.phecode.20231114.txt")
## expand by phecode
meddra.icd10  <- do.call(rbind, lapply(1:nrow(meddra.icd10), function(x) return(data.table(meddra.icd10[x,], phecode.single = strsplit(meddra.icd10$phecodes[x], "\\|")[[1]]))))
meddra.icd10  <- meddra.icd10[ phecode.single != ""]

## SNOMED-CT and ICD-10
snomed.icd10  <- fread("/sc-projects/sc-proj-computational-medicine/people/Maik/14_phecode_GWAS/19_phecode_disease_mapping/input/Mapping.SNOMEDCT.ICD10.phecode.20231114.txt")
## expand by phecode
snomed.icd10  <- do.call(rbind, lapply(1:nrow(snomed.icd10), function(x) return(data.table(snomed.icd10[x,], phecode.single = strsplit(snomed.icd10$phecodes[x], "\\|")[[1]]))))
snomed.icd10  <- snomed.icd10[ phecode.single != ""]

## MEDDRA and SNOMED
meddra.snomed <- fread("/sc-projects/sc-proj-computational-medicine/people/Maik/14_phecode_GWAS/19_phecode_disease_mapping/input/Mapping.MEDDRA.SNOMEDCT.phecode.20231114.txt") 
## expand by phecode
meddra.snomed <- do.call(rbind, lapply(1:nrow(meddra.snomed), function(x) return(data.table(meddra.snomed[x,], phecode.single = strsplit(meddra.snomed$phecodes[x], "\\|")[[1]]))))
meddra.snomed <- meddra.snomed[ phecode.single != ""]

## add phecodes to ADE mappings
ade.mapping   <- as.data.table(ade.mapping)
ade.mapping[, phecode.id := apply(ade.mapping[, c("notation", "source"), with=F], 1, function(x){
  
  ## get all mapping phecodes based on source
  if(x[2] == "SNOMEDCT"){
    ## get relevant mappings
    s1 <- unique(snomed.icd10[ snomed_ct == x[1]]$phecode.single)
    s2 <- unique(meddra.snomed[ snomed_ct == x[1]]$phecode.single)
  }else if(x[2] == "ICD10"){
    ## get relevant mappings
    s1 <- unique(snomed.icd10[ icd10 == x[1]]$phecode.single)
    s2 <- unique(meddra.icd10[ icd10 == x[1]]$phecode.single)
  }else{
    ## get relevant mappings
    s1 <- unique(meddra.icd10[ meddra == x[1]]$phecode.single)
    s2 <- unique(meddra.snomed[ meddra == x[1]]$phecode.single)
  }
  
  ## return results
  return(paste(sort(unique(c(s1,s2))), collapse = "|"))
  
})]
## split by phecode
ade.mapping.phecode <- ade.mapping[ phecode.id != ""]
ade.mapping.phecode <- lapply(1:nrow(ade.mapping.phecode), function(x){
  ## split by phecode
  return(data.table(ade.mapping.phecode[x,], phecode.id.single = strsplit(ade.mapping.phecode$phecode.id[x], "\\|")[[1]]))
})
## combine again
ade.mapping.phecode <- do.call(rbind, ade.mapping.phecode)
  
## add ADEs and indications to the list
drug.gene.QTL.rev[, reported_indication := sapply(id.drug, function(x) paste(sort(unique(ade.summary[compound_chembl_id == x & drug_effect_type == "indication"]$drug_effect)), collapse = "|"))]
drug.gene.QTL.rev[, reported_ade := sapply(id.drug, function(x) paste(sort(unique(ade.summary[compound_chembl_id == x & drug_effect_type != "indication"]$drug_effect)), collapse = "|"))]

## map indication to phecodes
drug.gene.QTL.rev[, phecode_indication_mapping := apply(drug.gene.QTL.rev[, c("beneficial.phecodes", "reported_indication")], 1, function(x){
  ## get all possible beneficial phecodes
  phecodes    <- strsplit(x[1], "\\|")[[1]]
  ## reported indications
  indications <- strsplit(x[2], "\\|")[[1]]
  print(phecodes)
  ## if any
  if(length(phecodes) > 0){
    ## map back to IDs
    phe.id <- sapply(phecodes, function(c) res.ldsc$id[which(res.ldsc$phenotype == c)])
    ## find possibly mapping ade
    ade.id <- unique(ade.mapping.phecode[ phecode.id.single %in% phe.id & ade %in% indications]$phecode.id.single)
    ## return possible matchings
    if(length(ade.id) > 0){
      return(paste(sapply(ade.id, function(c) res.ldsc$phenotype[which(res.ldsc$id == c)]), collapse = "|"))
    }else{
      return("")
    }
  }else{
    return("")
  }
})]

## map adverse reactions to phecodes
drug.gene.QTL.rev[, phecode_ade_mapping := apply(drug.gene.QTL.rev[, c("adverse.phecodes", "reported_ade")], 1, function(x){
  ## get all possible beneficial phecodes
  phecodes    <- strsplit(x[1], "\\|")[[1]]
  ## reported indications
  indications <- strsplit(x[2], "\\|")[[1]]
  print(phecodes)
  ## if any
  if(length(phecodes) > 0){
    ## map back to IDs
    phe.id <- sapply(phecodes, function(c) res.ldsc$id[which(res.ldsc$phenotype == c)])
    ## find possibly mapping ade
    ade.id <- unique(ade.mapping.phecode[ phecode.id.single %in% phe.id & ade %in% indications]$phecode.id.single)
    ## return possible matchings
    if(length(ade.id) > 0){
      return(paste(sapply(ade.id, function(c) res.ldsc$phenotype[which(res.ldsc$id == c)]), collapse = "|"))
    }else{
      return("")
    }
  }else{
    return("")
  }
})]

## check whether ade occurred among beneficial phecodes
drug.gene.QTL.rev[, phecode_benefical_ade_mapping := apply(drug.gene.QTL.rev[, c("beneficial.phecodes", "reported_ade")], 1, function(x){
  ## get all possible beneficial phecodes
  phecodes    <- strsplit(x[1], "\\|")[[1]]
  ## reported indications
  indications <- strsplit(x[2], "\\|")[[1]]
  print(phecodes)
  ## if any
  if(length(phecodes) > 0){
    ## map back to IDs
    phe.id <- sapply(phecodes, function(c) res.ldsc$id[which(res.ldsc$phenotype == c)])
    ## find possibly mapping ade
    ade.id <- unique(ade.mapping.phecode[ phecode.id.single %in% phe.id & ade %in% indications]$phecode.id.single)
    ## return possible matchings
    if(length(ade.id) > 0){
      return(paste(sapply(ade.id, function(c) res.ldsc$phenotype[which(res.ldsc$id == c)]), collapse = "|"))
    }else{
      return("")
    }
  }else{
    return("")
  }
})]

## column of ADE data from BNF was available
drug.gene.QTL.rev[, bnf_ade := ifelse(reported_ade != "", 1, 0)]

## based from a drugs perspective
write.table(drug.gene.QTL.rev, "Putative.drug.repurposing.adverse.effects.QTL.fused.20231116.txt", sep="\t", row.names=F)

## up-to-here --> import manually curated drug annotation table
drug.gene.QTL.rev <- as.data.table(readxl::read_excel("Putative.drug.repurposing.adverse.effects.QTL.fused.20231116.xlsx"))

## --> some examples for the manuscript <-- ##

## number of positive examples (some drugs have multiple targets)
length(unique(drug.gene.QTL.rev[ matching.indication != ""]$name)) ## n = 108
length(unique(drug.gene.QTL.rev[ matching.indication != ""]$name))/length(unique(drug.gene.QTL.rev$name))
## consistent effect directions
# length(unique(drug.gene.QTL.rev[ matching.indication != "" & effect_conflict == "no"]$name)) ## n = 72
## inconsistent effect direction
108 - length(unique(drug.gene.QTL.rev[ matching.indication != "" & effect_conflict == "yes" & tissue.diff == 0]$name)) ## n = 26
length(unique(drug.gene.QTL.rev[ matching.indication != "" & effect_conflict == "yes" & tissue.diff > 0]$id.drug)) ## n = 30

## positive examples with coherent effects
length(unique(drug.gene.QTL.rev[ !is.na(matching.indication) & effect_conflict == "no"]$name)) ## n = 63
## partially mixed effects
length(unique(drug.gene.QTL.rev[ !is.na(matching.indication) & effect_conflict == "yes"]$name)) ## n = 45
## confirmed ADEs
length(unique(drug.gene.QTL.rev[ !is.na(phecode_ade_mapping) ]$name)) ## n = 45



## SOAT1 example, 
exp(0.891795084816094); exp(0.891795084816094-1.96*0.146414249029377); exp(0.891795084816094+1.96*0.146414249029377)
exp(-3.30524517598905); exp(-3.30524517598905-1.96*0.542652677212506); exp(-3.30524517598905+1.96*0.542652677212506)

## how many drugs with at least one adverse effect, but no evidence for matching ADE
nrow(drug.gene.QTL.rev[ !is.na(adverse.phecodes) &  is.na(phecode_ade_mapping) & bnf_ade == 1 ])
## n = 50

## drugs with only beneficial effects
jj <- unique(drug.gene.QTL.rev[ !is.na(beneficial.phecodes) & is.na(adverse.phecodes)]$name)
## drop those that might have adverse effects via other routes
jj <- jj[-which(jj %in% unique(drug.gene.QTL.rev[ is.na(beneficial.phecodes) & !is.na(adverse.phecodes)]$name))]
## n = 305

## drugs with solely adverse effects
jj <- unique(drug.gene.QTL.rev[ is.na(beneficial.phecodes) & !is.na(adverse.phecodes)]$name)
## remove some with positive effects via other routes
jj <- jj[-which(jj %in% unique(drug.gene.QTL.rev[ !is.na(beneficial.phecodes) & is.na(adverse.phecodes)]$name))]
## n = 273

## drugs with both
length(unique(drug.gene.QTL.rev[ beneficial.phecodes != "" & adverse.phecodes != ""]$name)) ## n = 108

## most frequent 'side effects'
tail(sort(table(unlist(lapply(drug.gene.QTL.rev$adverse.phecodes, function(x) strsplit(x, "\\|")[[1]])))))
## Palpitations, Hypothyroidism,..

## how many high-confidence cis-eQTL/cis-pQTL assignments
res.all[, high_conf_cis := apply(res.all[, c( "eQTL_summary.coloc", "pQTL_summary.coloc", "candidate.gene.locus.r2"), with=F], 1, function(x){
  
  ## replace all missing values with "" to make things easier
  x[is.na(x)] <- ""
  
  ## go through and extract the list of possible target genes
  if(x[1] != ""){
    ## get all possible eGenes (prefer coloc over r2)
    e.gen  <- strsplit(x[1], "\\|\\|")[[1]]
    ## get all eGenes
    e.gen        <- as.data.frame(do.call(rbind, lapply(e.gen, function(k) strsplit(k, "\\(|\\)")[[1]][1:2])))
    ## add names
    names(e.gen) <- c("gene", "PP")
    ## convert PPs to numeric
    e.gen$PP     <- as.numeric(gsub("%", "", e.gen$PP))
    ## apply threshold
    e.gen$PP     <- subset(e.gen, PP >= .8)$gene
  }else{
    e.gen        <- NA
  }
  
  ## same for pQTL
  if(x[2] != ""){
    ## get all possible genes
    e.prot <- strsplit(x[2], "\\|\\|")[[1]]
    ## get all eGenes
    e.prot        <- as.data.frame(do.call(rbind, lapply(e.prot, function(k) strsplit(k, "\\(|\\)")[[1]][1:2])))
    ## add names
    names(e.prot) <- c("gene", "PP")
    ## convert PPs to numeric
    e.prot$PP     <- as.numeric(gsub("%", "", e.prot$PP))
    ## apply threshold
    e.prot$PP     <- subset(e.prot, PP >= .8)$gene
  }else{
    e.prot        <- NA
  }
  
  ## check for overlap with candidate gene assignment
  c.gen <- strsplit(x[3], "\\|")[[1]]
  
  ## return overlap
  return(paste(unique(intersect(c(e.gen, e.prot), c.gen)), collapse = "|"))
  
})]

## how many unique high-confidence assignments with QTL support
length(unique(res.all$high_conf_cis))
## n = 792 - 1


#########################################
####     positive examples drugs     ####
#########################################

## how many positive examples among drugs with QTL evidence
tmp         <- do.call(data.frame, aggregate(match.indication ~ id.drug, res.qtl.rev, function(x) sum(x)))
## how many indications are covered for those that also have a genetic finding (create intermediate data frame to ease mapping)
foo         <- lapply(1:nrow(res.ldsc), function(x){
  ## create output
  if(res.ldsc$drug.id[x] == ""){
    return(data.frame(id=res.ldsc$id[x], id.drug=""))
  }else{
    return(data.frame(id=res.ldsc$id[x], id.drug=strsplit(res.ldsc$drug.id[x], "\\|")[[1]]))
  }
}) 
foo         <- do.call(rbind, foo)
## reduce to phecodes with at least one genetic finding
foo         <- subset(foo, id %in% unique(res.all$id))
## indicate whether mapping could have been found in our data
tmp$phecode <- sapply(tmp$id.drug, function(x) paste(foo$id[which(foo$id.drug == x)], collapse = "|"))
## compress to unique brand names
tmp         <- merge(tmp, unique(res.qtl[, c("id.drug", "name")]))
## compress
tmp         <- lapply(unique(tmp$name), function(x){
  ## get the relevant data set
  jj <- subset(tmp, name == x)
  ## return
  if(nrow(jj) == 1){
    return(jj)
  }else{
    return(data.frame(id.drug = paste(jj$id.drug, collapse = "|"), match.indication=max(jj$match.indication), 
                      phecode=paste(unique(unlist(lapply(jj$phecode, function(c) strsplit(c, "\\|")[[1]]))), collapse = "|"), 
                      name=x))
  }
})
tmp         <- do.call(rbind, tmp)

## how many meds
nrow(tmp) ## n = 808 meds
## how many matching indications
nrow(subset(tmp, match.indication > 0)) ## n = 126
nrow(subset(tmp, match.indication == 0)) ## n = 559
## how many with possible evidence
nrow(subset(tmp, phecode != "" & match.indication == 0)) ## n = 508
nrow(subset(tmp, phecode == "")) ## n = 51

##########################################
####    adopt drug repurposing table  ####
##########################################

## adopt drug repurposing table, to have genetically validated drug effects next to 
## potentially new indications
drug.repurposing.updated <- drug.repurposing[, c("drug.id", "drug.locus.id", "name", "phenotype", "candidate.gene.locus.r2", "target.symbol", "type.gene", "drug.id.gwas", "drug.gwas.trial.matching", "drug.gwas.net.trial.matching", "match.indication", "blackBoxWarning", "description", "drugType", "hasBeenWithdrawn", "names.indication",
                                                 "ID", "ALLELE0", "ALLELE1", "A1FREQ", "BETA", "SE", "LOG10P", "related.to.index.condition")]
## split
tmp1                     <- drug.repurposing.updated[ match.indication == T]
tmp2                     <- drug.repurposing.updated[ match.indication == F]
## make unique (careful: some medications have different descriptions)
tmp2                     <- tmp2[ order(name, phenotype, drug.locus.id, description)]
tmp2[, ind := 1:.N, by = c("name", "phenotype", "drug.locus.id")]
tmp2                     <- tmp2[ ind == 1] ## n = 383
## create a new data set
drug.repurposing.updated <- merge(tmp2, unique(tmp1[, c("phenotype", "drug.locus.id", "drug.id.gwas", "drug.gwas.trial.matching", "drug.gwas.net.trial.matching", 
                                                        "ID", "ALLELE0", "ALLELE1", "A1FREQ", "BETA", "SE", "LOG10P")]), 
                                  by=c("ID", "drug.locus.id"), suffixes = c(".repurpose", ".indication"))

## keep only one indication for each 'new' drug
drug.repurposing.updated <- drug.repurposing.updated[ order(name, phenotype.repurpose, drug.locus.id, -LOG10P.indication)]
drug.repurposing.updated[, ind := 1:.N, by=c("name", "phenotype.repurpose", "drug.locus.id")]
## subset accordingly
drug.repurposing.updated <- drug.repurposing.updated[ ind == 1]

## write table to file
write.table(drug.repurposing.updated, "Updated.Drug.repurposing.table.20231201.txt", sep="\t", row.names=F)

