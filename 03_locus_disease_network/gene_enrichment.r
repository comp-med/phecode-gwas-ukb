######################################################
#### Perform gene enrichment on BiLouvain-        ####
#### Clustering results                           ####
####                                              ####
#### Hannah E. Schmidt                            ####
######################################################


rm(list=ls())
setwd("<path to file>")
options(stringsAsFactors = F)

## packages needed
install.packages("gprofiler2")
install.packages("fread")
install.packages('plyr')
require(data.table)
library(data.table)
## function to perform gene enrichment analysis
require('gprofiler2')



###### GENE LOCI ENRICHMENT #######

################################################
####           import relevant data         ####
################################################

## import LDSC Phecode GWAS file
res.ldsc  <- fread("<path to file>")

## import cluster file (result from BiLouvain Clustering), this file should not contain any nan in CoCluster_ID column
res.all    <- fread("<path to file>")

names(res.all)[names(res.all)=="V1"] <- "candidate.gene.locus.r2"
names(res.all)[names(res.all)=="V2"] <- "CoCluster ID"

clusters_more_than_three <- c()


#------------------------------------#
##--    run enrichment analysis   --##
#------------------------------------#

## test for enriched pathways for each CoCluster_ID (loci_cluster)
disease.pathways <- lapply(unique(res.all$CoCluster_ID), function(x){
  print("Cluster")
  print(x)
  ## get all relevant candidate genes
  tmp <- unique(res.all[ CoCluster_ID == x, candidate.gene.locus.r2])
  ## resolve ambiguous loci
  tmp <- unique(unlist(lapply(tmp, function(k) strsplit(k, "\\|\\|")[[1]])))
  if(length(tmp) >= 3){

  ## do the enrichment
  en.res <- gost(query = tmp, 
                 organism = "hsapiens", ordered_query = FALSE, 
                 multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                 measure_underrepresentation = FALSE, evcodes = TRUE, 
                 user_threshold = 0.05, correction_method = "gSCS", 
                 domain_scope = "annotated",
                 numeric_ns = "", sources = NULL, as_short_link = FALSE)
  ## return results
  return(data.table(Cluster=x, en.res$result))
  }
})

disease.pathways <- do.call(plyr::rbind.fill, disease.pathways)
## convert to data table
disease.pathways <- as.data.table(disease.pathways)
## drop findings with no enriched pathways

##drop intersection_size < 2
disease.pathways <- disease.pathways[disease.pathways$intersection_size>2]

##drop sort by p-values
disease.pathways <- disease.pathways[ !(is.na(p_value))]

#------------------------------------#
##--        write to file         --##
#------------------------------------#

## write relevant data to file to save enrichment analysis
write.table(disease.pathways[, -which(names(disease.pathways) == "parents"), with=F], "Results.Pathway.Enrichment.loci.txt", sep = "\t", row.names=F)




###### GENE PHENO ENRICHMENT #######

################################################
####           import relevant data         ####
################################################

## import LDSC Phecode GWAS file
res.ldsc  <- fread("<path to file>")
## import cluster file (result from BiLouvain Clustering)
res.all    <- fread("<path to file>")


clusters_more_than_three <- c()

#------------------------------------#
##--    run enrichment analysis   --##
#------------------------------------#

## test for enriched pathways for each phenotype_cluster
disease.pathways <- lapply(unique(res.all$phenotype_cluster), function(x){
  print("Cluster")
  print(x)
  ## get all relevant candidate genes
  tmp <- unique(res.all[phenotype_cluster == x, candidate.gene.locus.r2])
  ## resolve ambiguous loci
  tmp <- unique(unlist(lapply(tmp, function(k) strsplit(k, "\\|\\|")[[1]])))
  if(length(tmp) >= 3){
    
    ## do the enrichment
    en.res <- gost(query = tmp, 
                   organism = "hsapiens", ordered_query = FALSE, 
                   multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                   measure_underrepresentation = FALSE, evcodes = TRUE, 
                   user_threshold = 0.05, correction_method = "gSCS", 
                   domain_scope = "annotated",
                   numeric_ns = "", sources = NULL, as_short_link = FALSE)
    ## return results
    return(data.table(Cluster=x, en.res$result))
  }
})

disease.pathways <- do.call(plyr::rbind.fill, disease.pathways)
## convert to data table
disease.pathways <- as.data.table(disease.pathways)
## drop findings with no enriched pathways

##drop intersection_size < 2
disease.pathways <- disease.pathways[disease.pathways$intersection_size>2]

##drop sort by p-values
disease.pathways <- disease.pathways[ !(is.na(p_value))]

#------------------------------------#
##--        write to file         --##
#------------------------------------#

## write relevant data to file to save enrichment analysis
## write results to file
write.table(disease.pathways[, -which(names(disease.pathways) == "parents"), with=F], "Results.Pathway.Enrichment.pheno.txt", sep = "\t", row.names=F)
