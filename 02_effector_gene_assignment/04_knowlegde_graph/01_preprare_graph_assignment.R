###########################################################
#### use knowledge graph for candidate gene assignment ####
#### Maik Pietzner                                     ####
###########################################################

rm(list=ls())
setwd("<path to file>")
options(stringsAsFactors = F)
load(".RData")

## enable faster data reading
require(data.table)
## to work with graph objects
require(igraph)

###############################################
####      import phecode GWAS results      ####
###############################################

## import results
res.fine              <- fread("<path to file>")
## does already contain a column with closest genes ('closest.genes')

## import phecode summary
res.pheno             <- fread("<path to file>")

## import phecode EFO mapping
efo.phecodes          <- fread("<path to file>")
## add ICD-9 mappings
icd9.phecode          <- read.csv("<path to file>")
## add to phecode list
efo.phecodes$icd9     <- sapply(efo.phecodes$phecode, function(x){paste(subset(icd9.phecode, phecode == x)$icd9, collapse = "|")})
## careful: ICD10 does not have a decimal separator

## add american version of ICD-10 (CM)
icd10.cm.phecode      <- read.table("<path to file>", sep="\t", header=T)
## reduce to icd10
icd10.cm.phecode      <- subset(icd10.cm.phecode, flag == 10)
## add to phecode mapping
efo.phecodes$icd10.cm <- sapply(efo.phecodes$phecode, function(x){paste(subset(icd10.cm.phecode, phecode == x)$ICD, collapse = "|")})

###############################################
####             import primeKG            ####
###############################################

## import knowledge graph (https://www.biorxiv.org/content/10.1101/2022.05.01.489928v1)
prime.KG   <- fread("../download/kg.csv")
## import disease mappings
disease.KG <- fread("../download/disease_features.csv")
## same for durgs
drugs.KG   <- fread("../download/drug_features.csv")
## look at nodes
nodes.KG   <- fread("../download/nodes.csv")
## import edges
edges.KG   <- fread("../download/edges.csv")

## --> convert to graph <-- ##

## create graph
prime.KG.graph <- graph_from_data_frame(edges.KG[, c("x_index", "y_index", "relation", "display_relation")], directed = F, nodes.KG)

## --> extract subgraphs <-- ##

## import function to do so
source("../scripts/generate_subgraph.R")

## get node index
ii <- unique(subset(disease.KG, mondo_id %in% strsplit(efo.phecodes$mondo.id[which(efo.phecodes$phenotype == "Raynaud's syndrome")], "\\|")[[1]])$node_index)
## draw the network (look only at direct neighbours)
plot.pdf.sub(prime.KG.graph, c(ii, 2472)+1)
plot.pdf.sub(prime.KG.graph, 2472+1)

## some random networks
subset(nodes.KG, node_name == "INHBC")
plot.pdf.sub(prime.KG.graph, 58762+1)
subset(nodes.KG, node_name == "DKKL1")
plot.pdf.sub(prime.KG.graph, 3643+1)
subset(nodes.KG, node_name == "GRP")
plot.pdf.sub(prime.KG.graph, 11794+1)

###############################################
####       MONDO DB to map to phecodes     ####
###############################################

## function to import MONDO
mondo.db              <- ontologyIndex::get_OBO("../input/mondo.obo", extract_tags = "everything")
## get all xref entries into large mapping data frame

## do in parallel
require(doMC)
registerDoMC(12)

## loop through all mondo entries
mondo.map             <- mclapply(mondo.db$xref, function(x){
  
  ## separate source from entry
  x        <- lapply(x, function(c) strsplit(c, ":")[[1]])
  ## convert to data frame
  x        <- data.frame(do.call(rbind, x))
  ## report only if any
  if(nrow(x) > 0){
    ## add names
    names(x) <- c("source", "code")
    return(x)
  }
  
}, mc.cores = 12)
## combine
mondo.map             <- do.call(rbind, mondo.map)
## add MONDO code
mondo.map$mondo.id    <- gsub("\\.[0-9]*", "", row.names(mondo.map))
## create numeric identifier to map to primeKG
mondo.map$mondo.num   <- as.numeric(gsub("MONDO:", "", mondo.map$mondo.id))

## --> add to phecode mapping <-- ##

## go through each phecode
efo.phecodes$mondo.id <- sapply(1:nrow(efo.phecodes), function(x){
  
  ## get all possibly relevant identifiers
  icd9       <- strsplit(efo.phecodes$icd9[x], "\\|")[[1]]
  icd10.cm   <- strsplit(efo.phecodes$icd10.cm[x], "\\|")[[1]]
  snomed     <- strsplit(efo.phecodes$snomed[x], "\\|")[[1]]
  
  ## subset mondo list (ensure appropriate matching)
  mondo.i9   <- subset(mondo.map, source == "ICD9" & code %in% icd9)$mondo.num
  mondo.i10m <- subset(mondo.map, source == "ICD10CM" & code %in% icd10.cm)$mondo.num
  mondo.sno  <- subset(mondo.map, source == "SCTID" & code %in% snomed)$mondo.num
  
  ## return mapping
  return(paste(unique(c(mondo.i9, mondo.i10m, mondo.sno)), collapse = "|"))
  
})


## --> add MONDO ID to fine-mapping results <-- ##

## add to GWAS results
res.fine <- merge(res.fine, efo.phecodes[, c("phecode", "mondo.id")], by="phecode", all.x=T)
## careful not all map to a Mondo ID

###############################################
####  implement graph-gene-prioritization  ####
###############################################

## do in parallel
require(doMC)
registerDoMC(12)

## query each phecode to obtain candidate genes
tmp      <- mclapply(1:nrow(res.fine), function(x){
  

  ## proceed only if any
  if(mondo != ""){
    ## search whether there is a Mondo ID
    mondo   <- strsplit(res.fine$mondo.id[x], "\\|")[[1]]
    ## get node indices
    node.id <- unique(subset(disease.KG, mondo_id %in% mondo)$node_index)
    ## extract subgraph
    tmp.net <- ego(prime.KG.graph, order=2, nodes = node.id+1, mode = "all", mindist = 0)
    ## search for any genes
    tmp.net <- nodes.KG[ node_index %in% unlist(tmp.net) & node_type == "gene/protein"]
    ## test whether any of those is among the closest genes
    c.genes <- strsplit(res.fine$closest.genes[x], "\\|")[[1]]
    ## return overlap
    return(data.table(res.fine[x,], network.gene=paste(intersect(c.genes, tmp.net$node_name), collapse = "|")))
  }else{
    return(data.table(res.fine[x,], network.gene=""))
  }
  
}, mc.cores = 12)
## combine everything
tmp      <- do.call(rbind, tmp)

## rename and store information
res.fine <- tmp

## write to file
write.table(tmp, "<path to file>", sep="\t", row.names = F)

########################################################################
####                   Revision NatGen - 26/09/2023                 ####
########################################################################

#-----------------------------#
##-- check GO-term linkage --##
#-----------------------------#

## do in parallel
require(doMC)
registerDoMC(12)

## query each phecode to obtain candidate genes; update to include range of genes assigned
tmp      <- mclapply(1:nrow(res.fine), function(x){
  
  ## proceed only if any
  if(mondo != ""){
    ## search whether there is a Mondo ID
    mondo   <- strsplit(res.fine$mondo.id[x], "\\|")[[1]]
    ## get node indices
    node.id <- unique(subset(disease.KG, mondo_id %in% mondo)$node_index)
    ## extract subgraph
    tmp.net <- ego(prime.KG.graph, order=2, nodes = node.id+1, mode = "all", mindist = 0)
    ## search for any genes
    tmp.net <- nodes.KG[ node_index %in% unlist(tmp.net) & node_type == "gene/protein"]
    ## test whether any of those is among the closest genes
    c.genes <- strsplit(res.fine$closest.genes[x], "\\|")[[1]]
    ## return overlap
    return(data.table(res.fine[x,], network.gene=paste(intersect(c.genes, tmp.net$node_name), collapse = "|"), network.total=length(c.genes)))
  }else{
    return(data.table(res.fine[x,], network.gene=""))
  }
  
}, mc.cores = 12)
## combine everything
tmp      <- do.call(rbind, tmp)

## overall expansion
summary(tmp$network.total); sd(tmp$network.total)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00    4.00    9.00   13.23   18.00   74.00 

## overall after reducing to genes in close proximity
summary(sapply(tmp$network.gene, function(x) length(strsplit(x, "\\|")[[1]])))
sd(sapply(tmp$network.gene, function(x) length(strsplit(x, "\\|")[[1]])))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   0.000   0.000   1.537   2.000  32.000

