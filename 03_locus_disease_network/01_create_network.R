######################################################
#### Create a locus - phecode network             ####
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
require(igraph)

################################################
####           import relevant data         ####
################################################

## import collated results list
res.all        <- fread("<path to file>")

## drop MHC region
res.all        <- res.all[ group != "6_11"]

## import summary by phenotype
res.pheno      <- fread("<path to file>")

## get the colour code for phecode categories
cl.dat         <- aggregate(y.pos ~ category, res.pheno, function(x) c(min(x), max(x)))
cl.dat         <- do.call(data.frame, cl.dat)
names(cl.dat)  <- c("category", "start", "end")
cl.dat$mid     <- cl.dat$start + (cl.dat$end - cl.dat$start)*.5
## order
cl.dat         <- cl.dat[order(cl.dat$start),]
## add colour
cl.dat$cl      <- sapply(cl.dat$category, function(x) res.pheno$cl[which(res.pheno$category == x)][1])
## edit for colour labels as well
cl.dat$category[which(cl.dat$category == "")] <- "other"

################################################
####      create bipartite network data     ####
################################################

## --> packages needed <-- ##
require(igraph)

## create data frames to construct the network
e.tmp        <- res.all[, c("locus_id.r2", "phenotype", "BETA", "LOG10P")] 
## n=5956 edges

## now vertices
tmp1          <- data.frame(unique(res.all[, c("phenotype", "category")]), type="phenotype")
## rename
names(tmp1)   <- c("name", "category", "type")
tmp1$label    <- tmp1$name
## add locus data
tmp2          <- data.frame(locus_id.r2=unique(e.tmp$locus_id.r2), type="locus")
## add candidate gene annotation
tmp           <- aggregate(candidate.gene.locus.r2 ~ locus_id.r2,  res.all, function(x){
  paste(unique(unlist(lapply(x, function(k) strsplit(k, "\\|\\|")[[1]]))), collapse = "|")
})
## add gene annotation
tmp2          <- merge(tmp2, tmp)
## change names
names(tmp2)   <- c("name", "type", "label")
tmp2$category <- "locus"
## combine into one larger list
v.tmp         <- rbind(tmp1, tmp2)
## add colours
cl.dat        <- rbind(cl.dat, data.frame(category="locus", start=NA, end=NA, mid=NA, cl="white"))
## add to vertice table
v.tmp         <- merge(v.tmp, cl.dat[, c("category", "cl")])
## change order
v.tmp         <- v.tmp[, c("name", "type", "label", "cl", "category")]

## --> create the graph <-- ##
phecode.net   <- graph_from_data_frame(e.tmp, directed = F, v.tmp)

################################################
####        cluster bipartite graph         ####
################################################

## import bi-clustering results from Hannah
res.biclust                  <- fread("<path to file>")

## import corresponding enrichment results
enr.biclust                  <- fread("<path to file>")

## add column to biclust table with possible enriched pathways
res.biclust$enriched.pathway <- sapply(res.biclust$`CoCluster ID`, function(x){
  ## grep all associated pathways
  tmp <- enr.biclust[ Cluster == x]
  ## return if any
  if(nrow(tmp) > 0){
    ## order by p-value
    tmp <- tmp[ order(p_value)]
    return(paste(tmp$term_name, collapse = "|"))
  }else{
    return("")
  }
}) 

#---------------------------------#
##--  create induced subgraph  --##
#---------------------------------#

## get all nodes included in the biclustering
jj                  <- unlist(lapply(res.biclust$Phenotype, function(x) strsplit(x, ",")[[1]]))
## replace with names
jj                  <- sapply(jj, function(x) res.pheno$phenotype[which(res.pheno$id == x)])
## get the loci
ii                  <- unlist(lapply(res.biclust$Loci, function(x) strsplit(x, ",")[[1]]))

## create induced subgraph
bi.graph            <- induced_subgraph(phecode.net, c(ii, jj))

## add biclust community membership to each node (unfold)
tmp                 <- lapply(res.biclust$`CoCluster ID`, function(x){
  
  ## get respective row
  jj <- which(res.biclust$`CoCluster ID` == x)
  ## get phecodes
  ii <- strsplit(res.biclust$Phenotype[jj], ",")[[1]]
  ## rename
  ii <- sapply(ii, function(k) res.pheno$phenotype[which(res.pheno$id == k)])
  ## get loci
  kk <- strsplit(res.biclust$Loci[jj], ",")[[1]]
  ## return results
  return(data.table(cluster=x, name=c(ii, kk)))
  
})
## combine again
tmp                 <- do.call(rbind, tmp)

## add to the bigraph node attributes
V(bi.graph)$cluster <- sapply(V(bi.graph)$name, function(x) tmp$cluster[which(tmp$name == x)])

## export for plotting in cytoscape
edges.bilcust.cyto  <- as.data.table(igraph::as_data_frame(bi.graph))
## get locus to be always first
edges.bilcust.cyto[, locus := ifelse(edges.bilcust.cyto$from %in% res.all$locus_id.r2, edges.bilcust.cyto$from, edges.bilcust.cyto$to)]
edges.bilcust.cyto[, phenotype := ifelse(edges.bilcust.cyto$from %in% res.all$phenotype, edges.bilcust.cyto$from, edges.bilcust.cyto$to)]
## delete old estimates
edges.bilcust.cyto$BETA <- edges.bilcust.cyto$LOG10P <- NULL

#------------------------------------#
##-- add aligned effect estimates --##
#------------------------------------#

## define list of variants to be used (strongest for each r2 locus)
res.all             <- res.all[ order(locus_id.r2, -LOG10P)]
res.all[, ind := 1:.N, by="locus_id.r2"]
## get list of variants to be queried
tmp                 <- unique(res.all[ ind == 1, c("MarkerName", "ID", "locus_id.r2")])
## only those included in the network
tmp                 <- tmp[ locus_id.r2 %in% edges.bilcust.cyto$locus]
## get effect estimates
snp.effects         <- fread("<path to file>")
names(snp.effects)  <- c("misc", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA")
## create chromosome and phenotype
snp.effects$CHROM   <- as.numeric(gsub(".*:", "", snp.effects$misc)) 
snp.effects$id.phe  <- gsub("<path to file>/gwas_results/|\\.allchr.results.gz", "", gsub(":.*", "", snp.effects$misc))
## create MarkerName
snp.effects[, MarkerName := paste0("chr", as.numeric(gsub("X", 23, CHROM)), ":", as.numeric(GENPOS), "_", pmin(ALLELE0,ALLELE1), "_", pmax(ALLELE0,ALLELE1))]
## add locus
snp.effects         <- merge(tmp, snp.effects)
## add phenotype
snp.effects         <- merge(snp.effects, res.pheno[, c("id", "phenotype")], by.x="id.phe", by.y = "id")

## add to edges
edges.bilcust.cyto  <- merge(edges.bilcust.cyto, snp.effects[, c("locus_id.r2", "MarkerName", "ID", "BETA", "LOG10P", "phenotype")],
                             by.x=c("locus", "phenotype"), by.y=c("locus_id.r2", "phenotype"))

## get corresponding nodes
nodes.biclust.cyto  <- subset(v.tmp, name %in% c(edges.bilcust.cyto$locus, edges.bilcust.cyto$phenotype))
## n = 1127

## add cluster assignment
tmp                 <- lapply(res.biclust$`CoCluster ID`, function(x){
  
  ## get respective row
  jj <- which(res.biclust$`CoCluster ID` == x)
  ## get phecodes
  ii <- strsplit(res.biclust$Phenotype[jj], ",")[[1]]
  ## rename
  ii <- sapply(ii, function(k) res.pheno$phenotype[which(res.pheno$id == k)])
  ## get loci
  kk <- strsplit(res.biclust$Loci[jj], ",")[[1]]
  ## return results
  return(data.table(cluster=x, name=c(ii, kk)))
  
})
## combine again
tmp                 <- do.call(rbind, tmp)

## add to the nodes
nodes.biclust.cyto  <- merge(nodes.biclust.cyto, tmp)

#------------------------------------#
##--        write to file         --##
#------------------------------------#

## write relevant data to files to create cytoscape session
write.table(edges.bilcust.cyto, "Edges.biclustering.phecode.GWAS.20230208.txt", sep="\t", row.names = F, quote=F)
write.table(nodes.biclust.cyto, "Nodes.biclustering.phecode.GWAS.20230308.txt", sep="\t", row.names = F, quote=F)

################################################################################################################################
################################################################################################################################
########                                                    END OF SCRIPT                                               ######## 
################################################################################################################################
################################################################################################################################
