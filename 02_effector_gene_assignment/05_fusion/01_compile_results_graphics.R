######################################################
#### Compile summary of phecode GWAS discovery    ####
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

################################################
####       import different assignments     ####
################################################

## import overlap with OpenTargets
res.ot   <- fread("<path to file>")
## import overlap with GWAS catalog
res.gwas <- fread("<path to file>")
## import eQTL assignment (careful - misses MHC region)
res.eqtl <- fread("<path to file>")
## import pQTL assignment (careful - misses MHC region)
res.pqtl <- fread("<path to file>")
## import R2/CAFEH assignment
res.r2   <- fread("<path to file>")
## network assignment
res.net  <- fread("<path to file>")
## import VEP annotation
res.vep  <- fread("<path to file>")

################################################
####       combine into one large file      ####
################################################

#----------------------------------------#
##--       first two data sets        --##
#----------------------------------------#

## get common names to merge on
ii                 <- intersect(names(res.ot), names(res.gwas))
## drop some
ii                 <- ii[-which(ii %in% c("trait_reported", "study_id", "num_reported"))]

## merge both assignments
res.all            <- merge(res.ot, res.gwas, by=ii, suffixes = c(".ot", ".gwas"))

## generate novelty assignment
res.all[, unreported := ifelse(num_reported.gwas > 0, "GWAS_catalog", ifelse(num_reported.ot > 0, "OpenTargets", ifelse(group == "6_11", "unknown", "unreported")))] 

#----------------------------------------#
##--    add additional annotations    --##
#----------------------------------------#

## add eQTL results
res.all <- merge(res.all, res.eqtl[, c("id", "MarkerName", "eQTL_summary.coloc", "eQTL_summary.r2")], by=c("id", "MarkerName"), all.x=T)
## add pQTL results
res.all <- merge(res.all, res.pqtl[, c("id", "MarkerName", "pQTL_summary.coloc", "pQTL_summary.r2")], by=c("id", "MarkerName"), all.x=T)
## add r2 linkage and key CAFEH output
res.all <- merge(res.all, res.r2[, c("id", "MarkerName", "R2.group", "locus_id.r2", "top_component", "p_active.top.susie", "locus_id.cafeh")], by=c("id", "MarkerName"), all.x=T)

#----------------------------------------#
##--    add functional annotation     --##
#----------------------------------------#

## define names to be kept
jj             <- c("lead_MarkerName", "id", "lead_CHROM", "lead_GENPOS", "lead_ID", "lead_cs", "lead_Consequence", "lead_IMPACT", "lead_SYMBOL", "lead_Gene", "lead_CADD_PHRED",
                    "proxy_R2", "proxy_MarkerName", "proxy_GENPOS", "proxy_ID", "proxy_pip", "proxy_Consequence", "proxy_IMPACT", "proxy_Gene", "proxy_SYMBOL", "proxy_CADD_PHRED")
## reduce
tmp            <- res.vep[, ..jj]
res.all[ group != "6_11" & !(MarkerName %in% tmp$lead_MarkerName)]
## keep only those included in res.all and combine again afterwards
tmp            <- tmp[lead_MarkerName %in% res.all$MarkerName]

## --> get variants in the MHC region and those missing <-- ##

## get information on variants in MHC region as well
all.vep        <- fread("<path to file>")
## create MarkerName
all.vep[, MarkerName:=paste0("chr", as.numeric(gsub("X", 23, CHROM)), ":", as.numeric(POS), "_", pmin(REF,ALT), "_", pmax(REF,ALT))]
## subset to variants in MHC region (ID not good to merge on, create MarkerName instead)
all.vep        <- all.vep[ MarkerName %in% unique(res.all[ !(MarkerName %in% tmp$lead_MarkerName), MarkerName])]
## get only what is really needed
all.vep        <- all.vep[, c("MarkerName", "CHROM", "POS", "Consequence", "IMPACT", "Gene", "SYMBOL", "CADD_PHRED"), with=F]
## change names to merge with other assignments
names(all.vep) <- paste0("lead_", names(all.vep))
## add phecode IDs
tmp2           <- merge(all.vep, res.all[ !(MarkerName %in% tmp$lead_MarkerName), c("MarkerName", "id")], by.x="lead_MarkerName", by.y="MarkerName")
tmp2           <- as.data.table(tmp2)

## --> combine <-- ##

## header to be kept
jj <- c("lead_Consequence", "lead_IMPACT", "lead_SYMBOL", "lead_Gene", "lead_CADD_PHRED",
        "proxy_R2", "proxy_MarkerName", "proxy_GENPOS", "proxy_ID", "proxy_pip", "proxy_Consequence",
        "proxy_IMPACT", "proxy_Gene", "proxy_SYMBOL", "proxy_CADD_PHRED")

## do in parallel
require(doMC)
registerDoMC(12)
## run as loop
foo            <- mclapply(1:nrow(res.all), function(x){
  
  ## get lead MarkerName and phecode ID
  phe <- res.all$id[x]
  mn  <- res.all$MarkerName[x]
  
  ## search for matching pattern in vep output
  vep <- res.vep[lead_MarkerName == mn & id == phe]
  
  ## add depending on matching
  if(nrow(vep) == 1){
    ## take only proxies forward with r2_proxy > .6
    if(vep$proxy_R2 >= .6){
      return(data.table(res.all[x, ], vep[, ..jj]))
    }else{
      return(data.table(res.all[x, ], vep[, c("lead_Consequence", "lead_IMPACT", "lead_SYMBOL", "lead_Gene", "lead_CADD_PHRED"), with=F]))
    }
    
  }else{
    
    ## ckeck whether within MHC region
    if( res.all$group[x] == "6_11" ){
      ## add minimum for lead variant only
      return(data.table(res.all[x,], tmp2[ id == phe & lead_MarkerName == mn, c("lead_Consequence", "lead_IMPACT", "lead_SYMBOL", "lead_Gene", "lead_CADD_PHRED"), with=F]))
    }else{
      ## create file name
      fname <- paste0("credible.set.", gsub("\\|", ".", res.all$locus_id[x]), ".txt")
      ## fuse results for proxy signals with ambiguous lead signals
      return(data.table(res.all[x,], tmp2[ id == phe & lead_MarkerName == mn, c("lead_Consequence", "lead_IMPACT", "lead_SYMBOL", "lead_Gene", "lead_CADD_PHRED"), with=F], 
                        res.vep[ lead_cs == res.all$cs[x] & file_name == fname, c("proxy_R2", "proxy_MarkerName", "proxy_GENPOS", "proxy_ID", "proxy_pip", "proxy_Consequence",
                                                                                  "proxy_IMPACT", "proxy_Gene", "proxy_SYMBOL", "proxy_CADD_PHRED")]))
    }
  }
  
}, mc.cores = 12) 
## combine everything
foo            <- do.call(plyr::rbind.fill, foo)
## convert to data table
foo            <- as.data.table(foo)

#----------------------------------------#
##--           add closest gene       --##
#----------------------------------------#

## add closest gene
require(biomaRt)

## get data on build 37
gene.ensembl                       <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37) 

## obtain a list of all protein coding genes
tmp.genes                          <- getBM(attributes = c('chromosome_name', 'start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype', 'transcription_start_site'),
                                            filters = c('chromosome_name'),
                                            values = list(c(1:22, "X")),
                                            mart = gene.ensembl)
## convert to data table to ease query
tmp.genes                          <- as.data.table(tmp.genes)
## convert some entries to numeric
tmp.genes$start_position           <- as.numeric(tmp.genes$start_position)
tmp.genes$end_position             <- as.numeric(tmp.genes$end_position)
tmp.genes$transcription_start_site <- as.numeric(tmp.genes$transcription_start_site)

## loop through each entry and define closest gene (protein coding)
tmp          <- lapply(1:nrow(foo), function(x){
  
  print(x)
  
  ## get the gene of potential interest (2MB window)
  tmp           <- tmp.genes[ chromosome_name == ifelse(foo$CHROM[x] == 23, "X", foo$CHROM[x]) & start_position >= foo$GENPOS[x]-1e6 & end_position <= foo$GENPOS[x]+1e6 ]
  
  ## test whether enough protein coding gens are in the region
  if(nrow(subset(tmp, gene_biotype %in% c("protein_coding", "processed_transcript"))) > 0){
    ## restrict to protein encoding genes for now
    tmp <- subset(tmp, gene_biotype %in% c("protein_coding", "processed_transcript"))
  }
  
  ## compute distance to any gene (gene body not TSS)
  tmp$dist.body <- apply(tmp[, c("start_position", "end_position")], 1, function(k) min(abs(k-foo$GENPOS[x])))
  ## compute distance to TSS
  tmp$dist.tss  <- abs(tmp$transcription_start_site - foo$GENPOS[x])
  ## sort by distance to gene body
  tmp           <- tmp[order(tmp$dist.body), ]
  ## get position for TSS
  jj            <- which.min(tmp$dist.tss)
  
  ## return three different entries
  return(data.table(foo[x, ], 
                    closest.gene.body=paste0(tmp$external_gene_name[1], " (", tmp$dist.body[1], ") - ", tmp$ensembl_gene_id[1]),
                    closest.gene.tss=paste0(tmp$external_gene_name[jj], " (", tmp$dist.tss[jj], ") - ", tmp$ensembl_gene_id[jj]),
                    closest.genes=paste(unique(subset(tmp, dist.body <= 5e5)$external_gene_name), collapse = "|")))
  
  
  
})
## combine everything
tmp     <- do.call(rbind, tmp)
save.image()

## re-assign
res.all <- tmp

#----------------------------------------#
##--         add network gene         --##
#----------------------------------------#

## add network gene assignment
res.all <- merge(res.all, res.net[, c("id", "MarkerName", "network.gene")])

################################################
####           gene prioritization          ####
################################################

## get description of VEP terms
rank.cons               <- data.frame(consequence=c("transcript_ablation", "splice_acceptor_variant", "splice_donor_variant", "stop_gained", "frameshift_variant", "stop_lost", "start_lost", "transcript_amplification", "inframe_insertion", "inframe_deletion",
                                                    "missense_variant", "protein_altering_variant", "splice_region_variant", "incomplete_terminal_codon_variant", "start_retained_variant", "stop_retained_variant",
                                                    "synonymous_variant", "coding_sequence_variant", "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_exon_variant", "intron_variant",
                                                    "NMD_transcript_variant", "non_coding_transcript_variant", "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant",
                                                    "regulatory_region_ablation", "regulatory_region_amplification", "feature_elongation", "regulatory_region_variant", "feature_truncation"), 
                                      rank=1:35)
## rank ≤ 12 considered functional

## create a variable for a possible functional gene
res.all$functional.gene <- apply(res.all[, c("lead_Consequence", "lead_SYMBOL", "proxy_Consequence", "proxy_SYMBOL"), with=F], 1, function(x){
  
  ## get possible consequence for both
  lead  <- strsplit(x[1], "&")[[1]]
  proxy <- strsplit(x[3], "&")[[1]]
  ## replace with numeric entry; add 36 to avoid problems with NA
  lead  <- c(36, rank.cons$rank[rank.cons$consequence %in% lead])
  proxy <- c(36, rank.cons$rank[rank.cons$consequence %in% proxy])
  ## get the maximum
  lead  <- min(lead, na.rm = T)
  proxy <- min(proxy, na.rm = T)
  ## report back
  if(lead <= 12 | proxy <= 12){
    if(lead <= proxy){
      return(x[2])
    }else{
      return(x[4])
    }
  }else{
    return("")
  }
})

#----------------------------------------#
##--     temporary gene assignment    --##
#----------------------------------------#

## go through each entry separately first
res.all$candidate.gene <- apply(res.all[, c("eQTL_summary.coloc", "pQTL_summary.coloc", "v2g.all", "closest.gene.body", "closest.gene.tss", "network.gene", "functional.gene")], 1, function(x){

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
    e.gen$PP     <- ifelse(e.gen$PP > .8, 1, ifelse(e.gen$PP > .5, .5, 0))
  }else{
    e.gen        <- data.frame(gene="", PP=NA)
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
    e.prot$PP     <- ifelse(e.prot$PP > .8, 1, ifelse(e.prot$PP > .5, .5, 0))
  }else{
    e.prot        <- data.frame(gene="", PP=NA)
  }
  
  ## OpenTarget annotation
  if(x[3] != ""){
    ## get all possible genes
    ot.gen        <- strsplit(x[3], "\\|")[[1]]
    ## get all eGenes
    ot.gen        <- as.data.frame(do.call(rbind, lapply(ot.gen, function(k) strsplit(k, "_")[[1]][c(1,3)])))
    ## add names
    names(ot.gen) <- c("gene", "PP")
    ## convert PPs to numeric
    ot.gen$PP     <- as.numeric(ot.gen$PP)
    ## apply threshold
    ot.gen$PP     <- ifelse(ot.gen$PP > .25, 1, ifelse(ot.gen$PP > .1, .5, 0))
  }else{
    ot.gen        <- data.frame(gene="", PP=NA)
  }
  
  ## same for network
  if(x[6] != ""){
    net.gen <- strsplit(x[6], "\\|")[[1]]
    net.gen <- data.frame(gene=net.gen, PP = 1/length(net.gen))
  }else{
    net.gen <- data.frame(gene="", PP=NA)
  }
  
  ## combine all three sets
  cand.gen    <- merge(e.gen, e.prot, by="gene", all=T, suffixes = c(".eQTL", ".pQTL"))
  cand.gen    <- merge(cand.gen, ot.gen, by="gene", all=T)
  ## network gene
  cand.gen    <- merge(cand.gen, net.gen, by="gene", suffixes = c(".ot", ".net"), all=T)
  ## add closest gene (body)
  cand.gen    <- merge(cand.gen, data.frame(gene=strsplit(x[4], " ")[[1]][1], PP=.5), by="gene", all=T)
  cand.gen    <- merge(cand.gen, data.frame(gene=strsplit(x[5], " ")[[1]][1], PP=.5), by="gene", all=T, suffixes = c(".body", ".tss"))
  ## add functional gene
  cand.gen    <- merge(cand.gen, data.frame(gene=x[7], PP=2), all=T, by="gene")
  ## delete possible missing genes
  cand.gen    <- subset(cand.gen, gene != "")
  
  ## compute final score
  cand.gen$sc <- apply(cand.gen[,-1], 1, sum, na.rm=T)
  ## get all meeting this score
  cand.gen    <- subset(cand.gen, sc == max(cand.gen$sc, na.rm=T))
  
  ## return top genes
  return(paste(apply(cand.gen[, c("gene", "sc")], 1, function(c) paste0(c[1], " (", c[2], ")")), collapse = "||"))
 
})

## divide gene and score
res.all$candidate.gene.score <- as.numeric(gsub(".*\\(|\\).*", "", res.all$candidate.gene))
res.all$candidate.gene       <- gsub("\\s*\\([^\\)]+\\)", "", res.all$candidate.gene)

## compute consensus gene by locus_id
tmp                          <- lapply(unique(res.all$locus_id.r2), function(x){
  
  ## treat MHC region differently
  if(is.na(x)){
    return(data.table(res.all[is.na(locus_id.r2), ], candidate.gene.locus.r2=res.all[is.na(locus_id.r2), candidate.gene]))
  }else{
    ## get all the relevant signals
    tmp <- res.all[ locus_id.r2 == x ]
    ## get the strongest candidates
    jj  <- tmp[ candidate.gene.score == max(tmp$candidate.gene.score), candidate.gene]
    
    ## return compressed list
    return(data.table(tmp, candidate.gene.locus.r2=paste(sort(unique(unlist(lapply(jj, function(k) strsplit(k, "\\|\\|")[[1]])))), collapse = "|")))
  }
  
})
## combine everything
tmp                         <- do.call(rbind, tmp)
## do assignment
res.all                     <- tmp

################################################
####              overall summary           ####
################################################

## import results from LD-score regression
res.ldsc         <- fread("<path to file>")
## compute overall number of participants
res.ldsc$n.total <- res.ldsc$cases + res.ldsc$controls

#----------------------------------------#
##-- add information on previous GWAS --##
#----------------------------------------#

## load mapping from 
gwas.mapping        <- fread("<path to file>")

## add to LDSC results
res.ldsc            <- merge(res.ldsc, gwas.mapping[, c("id", "m.freq", "icd10", "snomed", "ctv3", "read2", "GWAS.cat", "efo", "efo.description", "GWAS.collated")])

## import information on phecode resource composition
phecode.resource <- fread("<path to file>")

## add selected information to res.ldsc
res.ldsc         <- merge(res.ldsc, phecode.resource[, c("phecode", "primary_care", "all", "wo.primary_care", "case.ratio.pcare")])

## add case ratio also to SNV results
res.all          <- merge(res.all, phecode.resource[, c("phecode", "case.ratio.pcare")])
## create new variable
res.all[, pcare.most := ifelse(case.ratio.pcare >= .5, "yes", "no")]

################################################
####    create refined novelty assignment   ####
################################################

## map previous GWAS to results file and test whether reporting matched
tmp <- lapply(1:nrow(res.all), function(x){
  
  ## get the relevant phecode
  phe <- res.all$phecode[x]
  ## get all mapping GWAS studies, if any
  phe <- strsplit(res.ldsc$GWAS.collated[which(res.ldsc$phecode == phe)], "\\|")[[1]]
  ## proceed only if any
  if(length(phe) > 0){
    ## get all GCST codes
    phe <- unlist(lapply(phe, function(k) strsplit(k, " ")[[1]][1]))
    ## get all mapping GWAS entries
    gws <- strsplit(res.all$study_id.gwas[x], "\\|\\|")[[1]]
    ## return what might be mapping
    return(data.table(res.all[x,], unreported.refined=paste(intersect(phe, gws), collapse = "|")))
  }else{
    ## return with no added entry
    return(data.table(res.all[x,], unreported.refined=""))
  }
  
  
})
## combine again
tmp <- do.call(rbind, tmp)

## create additional column
tmp$unreported.trait.mapping <- tmp$unreported.refined
tmp$unreported.refined       <- ifelse(tmp$unreported.refined == "", 0, 1)

## reassign
res.all                      <- tmp

################################################
####       derive loci characteristics      ####
################################################

## create overlap plot for:
## -- novelty
## -- primary / secondary signal
## -- gwas-significance 
## -- colour by phecode categories
## -- allele frequency

## --> create column for lead or secondary signal in the region <-- ##
res.all <- as.data.table(res.all)
res.all <- res.all[order(locus_id, -LOG10P)]
## add indicator
res.all[, signal_locus := 1:.N, by="locus_id"]
## indicator whether lead or secondary signal
res.all[, sentinel := ifelse(signal_locus == 1, 1, 0)]

## --> create column for gwas-significance <-- ##
res.all$gw.sig.corrected <- ifelse(res.all$LOG10P > -log10(5e-8/nrow(res.ldsc)), 1, 0)

## --> indicator for allele frequency <-- ##
res.all$maf.category     <- ifelse(res.all$A1FREQ >= .99 | res.all$A1FREQ <= .01, 1, 0)

## --> create simplified columns for novelty <-- ##
res.all$reported.gwas.cat.mapping <- ifelse(res.all$unreported == "GWAS_catalog" & res.all$unreported.refined == 1, 1, 0) 
res.all$reported.gwas.cat.diff    <- ifelse(res.all$unreported == "GWAS_catalog" & res.all$unreported.refined == 0, 1, 0) 
res.all$reported.open.tar         <- ifelse(res.all$unreported == "OpenTargets", 1, 0)

#-------------------------------#
##--      create the plot    --##
#-------------------------------#

## edit one group variable
res.all$category[which(res.all$category == "")] <- "other"
## generate unique MarkerName - phecode variable
res.all$snp.phecode.id                          <- paste(res.all$MarkerName, res.all$id, sep="$")

## load function to do so
source("../scripts/plot_overlap.R")

## third version by whether condition is based mostly on primary care
pdf("../graphics/Overlap.SNP.Phecode.associations.primary.care.20240406.pdf", width = 3.15, height = 3.15)
# png("../graphics/Overlap.SNP.Phecode.associations.primary.care.20221129.png", width = 8, height = 8, units = "cm", res=900)
par(mar=c(.1,4,.5,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, xaxs="i", yaxs="i")
## define layout
layout(matrix(1:2, 2, 1), heights = c(.6,.4))

## try
plot.overlap(res.all, "snp.phecode.id", 
             c("gw.sig.corrected", "sentinel", "maf.category", "reported.gwas.cat.mapping", "reported.gwas.cat.diff", "reported.open.tar"),
             c("corrected GW-sig.", "Regional sentinel", "MAF<1%", "Replication\nGWAS catalog", "New trait\nGWAS catalog", "Reported\nOpen Targets"),
             "pcare.most", data.frame(pcare.most=c("no", "yes"), cl=c(colorspace::lighten("#F95700", .5), "#F95700")))

dev.off()


## write to file
write.table(res.all, "<path to file>", sep="\t", row.names=F)

#---------------------------------#
##--  import credible set size --##
#---------------------------------#

## do in parallel
require(doMC)
registerDoMC(12)

## import size of credible sets
tmp     <- mclapply(1:nrow(res.all), function(x){
  
  ## omit MHC region
  if(res.all$group[x] != "6_11"){
    
    ## get what is needed
    jj <- paste("credible.set", res.all$id[x], res.all$CHROM[x], res.all$region_start_collated[x], res.all$region_end_collated[x], "txt",sep=".")
    ## import the data
    jj <- fread(paste0("<path to file>", jj))
    ## return size of the relevant credible set
    return(data.table(res.all[x, ], cs.size = nrow(jj[ cs == res.all$cs[x]])))
  }else{
    return(data.table(res.all[x, ], cs.size=NA))
  }
}, mc.cores = 12)
## combine again
tmp     <- do.call(rbind, tmp)
## reassign
res.all <- tmp 

################################################
####              MAF vs beta plot          ####
################################################

## colour vector for novelty
cl.rep         <- data.frame(unreported=c("GWAS_catalog", "OpenTargets", "unreported", "unknown"),
                             cl.rep=c("grey80", "grey60", "grey10", "grey70"))

## add colour for novelty
res.all        <- merge(res.all, cl.rep, by="unreported")

## the actual plot
pdf("<path to file>", width = 3.15, height = 3.15)
par(mar=c(1.5,1.5,1,.5), mgp=c(.6,0,0), cex.axis=.5, cex.lab=.5, tck=.01, lwd=.5)

## empty plot
plot(c(0,.5), range(abs(res.all$BETA)), type="n", xlab="Minor allele frequency",
     ylab="Absolute effect size per minor allele [log(odds)]", xaxt="n", yaxt="n")
## add axis
axis(1, lwd=.5); axis(2, lwd=.5)
## add points for reported signals
tmp <- res.all[ unreported != "unreported"]
points(ifelse(tmp$A1FREQ > .5, 1-tmp$A1FREQ, tmp$A1FREQ), abs(tmp$BETA), pch=20,
       col=tmp$cl.rep, cex=.3)
## add unreported ones
tmp <- res.all[ unreported == "unreported"]
points(ifelse(tmp$A1FREQ > .5, 1-tmp$A1FREQ, tmp$A1FREQ), abs(tmp$BETA), pch=20,
       col=tmp$cl.rep, cex=.3)

## highlight some examples
tmp <- res.all[ abs(res.all$BETA) > 4]
## order
tmp <- tmp[order(-abs(BETA))]
## get plotting parameters
pm  <- par("usr")
## add labels
text(.05, seq(pm[4]-(pm[4]-pm[3])*.02, pm[4]-(pm[4]-pm[3])*.5, length.out=nrow(tmp)), 
     labels=paste(tmp$ID, tmp$candidate.gene, tmp$phenotype, sep=" - "),
     pos=4, offset=0, cex=.35)
## add arrows
arrows(.05, seq(pm[4]-(pm[4]-pm[3])*.02, pm[4]-(pm[4]-pm[3])*.5, length.out=nrow(tmp)),
       .025, seq(pm[4]-(pm[4]-pm[3])*.02, pm[4]-(pm[4]-pm[3])*.5, length.out=nrow(tmp)),
       length = 0, lwd=.3)
arrows(.025, seq(pm[4]-(pm[4]-pm[3])*.02, pm[4]-(pm[4]-pm[3])*.5, length.out=nrow(tmp)),
       ifelse(tmp$A1FREQ > .5, 1-tmp$A1FREQ, tmp$A1FREQ), abs(tmp$BETA),
       length = 0, lwd=.3)

## add legend
legend(pm[1], pm[4]+(pm[4]-pm[3])*.075, pch=20, pt.cex=.8, lty=0, bty="n", cex=.5,
       legend=cl.rep$unreported, col=cl.rep$cl.rep, xpd=NA, ncol=4)

dev.off()

################################################
####          plot overall results          ####
################################################

## generate genetic coordinates
gen.coord      <- data.table::fread("zcat <path to file>/gwas_results/date_1000.allchr.results.gz", sep=" ", header=T, select=c("CHROM", "GENPOS", "ID"))
## sort and add ordering
gen.coord      <- gen.coord[order(gen.coord$CHROM, gen.coord$GENPOS),]
gen.coord[, x.pos:=1:.N,]

## add to the fine-mapping results
res.all        <- merge(res.all, gen.coord[, c("ID", "x.pos")], by="ID")

## generate data for Chromosomes
chr.dat        <- merge(aggregate(x.pos ~ CHROM, gen.coord, min), aggregate(x.pos ~ CHROM, gen.coord, max), suffixes = c(".min", ".max"), by="CHROM") 
names(chr.dat) <- c("chr", "start", "end")
chr.dat$mid    <- chr.dat$start + (chr.dat$end - chr.dat$start)*.5

## import summary by phenotype
res.pheno      <- fread("<path to file>")

## colour vector for novelty
cl.rep         <- data.frame(unreported=c("GWAS_catalog", "OpenTargets", "unreported", "unknown"),
                             cl.rep=c("grey80", "grey60", "grey10", "grey70"))

## add colour for novelty
res.all        <- merge(res.all, cl.rep, by="unreported")
## add y-category
res.all        <- merge(res.all, res.pheno[, c("id", "y.pos")], by="id")

## add information whether or not a GWAS exists to SNP results
res.all        <- merge(res.all, res.ldsc[, c("id", "gwas")], by="id")
res.all        <- as.data.table(res.all)

## get new phecode ordering based on whether or not a GWAS existed
tmp              <- res.ldsc[, c("id", "gwas", "phecode", "category")]
tmp              <- tmp[order(gwas, phecode)]
tmp$y.pos.v2     <- 1:nrow(tmp)

## add new order
res.all          <- merge(res.all, tmp)

## create new set to colour phecode categories
cl.dat2          <- aggregate(y.pos.v2 ~ gwas + category, tmp, function(x) c(min(x), max(x)))
cl.dat2          <- do.call(data.frame, cl.dat2)
names(cl.dat2)   <- c("gwas", "category", "start", "end")
cl.dat2$mid      <- cl.dat2$start + (cl.dat2$end - cl.dat2$start)*.5
## order
cl.dat2          <- cl.dat2[order(cl.dat2$start),]
## add colour
cl.dat2$cl       <- sapply(cl.dat2$category, function(x) res.ldsc$cl[which(res.ldsc$category == x)][1])

#----------------------------------------#
##--         2D Manhattan plot        --##
#----------------------------------------#

## indicator, whether mostly diagnosed in primary care
phecode.resource[, pcare.most := ifelse(case.ratio.pcare >= .5, "yes", "no")]

## get new phecode ordering based on whether or not a GWAS existed
tmp              <- phecode.resource[, c("id", "pcare.most", "phecode", "category")]
tmp              <- tmp[order(pcare.most, phecode)]
tmp$y.pos.v2     <- 1:nrow(tmp)

## add new order
res.all          <- merge(res.all, tmp[, c("id", "y.pos.v2")], by="id")

## create new set to colour phecode categories
cl.dat2          <- aggregate(y.pos.v2 ~ pcare.most + category, tmp, function(x) c(min(x), max(x)))
cl.dat2          <- do.call(data.frame, cl.dat2)
names(cl.dat2)   <- c("pcare.most", "category", "start", "end")
cl.dat2$mid      <- cl.dat2$start + (cl.dat2$end - cl.dat2$start)*.5
## order
cl.dat2          <- cl.dat2[order(cl.dat2$start),]
## add colour
cl.dat2$cl       <- sapply(cl.dat2$category, function(x) res.ldsc$cl[which(res.ldsc$category == x)][1])

#----------------------------------------#
##--         2D Manhattan plot        --##
#----------------------------------------#

## collapse data by region instead of r2 locus
region.count <- lapply(unique(res.all$group), function(x){
  
  ## get the relevant list
  tmp <- res.all[ group == x]
  ## return relevant information
  return(data.table(group=x, n.assoc=length(unique(tmp$id)), 
                    ## prioritize IL6R
                    # top.gene=ifelse("IL6R" %in% tmp$candidate.gene.locus.r2, "IL6R", names(sort(table(tmp$candidate.gene.locus.r2), decreasing = T))[1]),
                    top.gene=names(sort(table(tmp$candidate.gene.locus.r2), decreasing = T)[1]),
                    n.categories=length(unique(tmp$category)), x.pos=mean(tmp$x.pos)))
  
})
## combine
region.count <- do.call(rbind, region.count)

## graphic
pdf("<path to file>", width = 7.1, height = 5.1)
par(mar=c(.1,1.5,3.5,.5), mgp=c(.4,0,0), cex.axis=.5, cex.lab=.5, tck=-.01, lwd=.5, xaxs="i", yaxs="i")
## define layout
layout(matrix(1:2, 2, 1), heights = c(.3,.7))

#---------------------------------#
##-- upper part for pleiotropy --##
#---------------------------------#

## empty plot
plot(c(chr.dat$start[1], chr.dat$end[23]), c(0,50), xaxt="n", yaxt="n",
     xlab="", ylab="Number\nassociated diseases", type="n")
## add axis
axis(2, lwd=.5)
## rectanlges to divide chromosomes
rect(chr.dat$start, 0, chr.dat$end, 50, col=c("white", "grey90"), border = NA)

## add regional estimates
arrows(region.count$x.pos, 0, region.count$x.pos, ifelse(region.count$n.assoc >= 50, 50, region.count$n.assoc),  length = 0, lwd=.5, col="grey50")

## add legend
pm <- par("usr")
legend(pm[1]+(pm[2]-pm[1])*.23, pm[4], ncol=1, bty="n", lty=0, pch=21, pt.bg="white", cex=.4, pt.lwd=.3,
       pt.cex=rev(c(1,5,10,15))/15, legend = rev(c(1,5,10,15)), title = "ICD-10 chapters")

## add genes with at least 5 annotations
tmp          <- subset(region.count, n.assoc >= 12)
## edit one name
tmp$top.gene <- gsub("BPNT1|RP11-95P13.2|SLC30A10|ZC3H11B", "BPNT1*", tmp$top.gene, fixed = T)
tmp$top.gene <- ifelse(tmp$n.assoc >= 50, paste0(tmp$top.gene, " (", tmp$n.assoc, ")"), tmp$top.gene)
tmp          <- tmp[order(tmp$x.pos), ]
## define plotting coordinates
pm  <- par("usr")

## define the coordinates of plotting the boxes
p.ii <- tmp$x.pos
## now get the distance between each point
d.ii <- p.ii[-1] - p.ii[-length(p.ii)]
## width of one box
w    <- strwidth("B")*1.7
## get the gaps
g.ii <- c(0,which(d.ii > w), length(p.ii))
for(j in 2:length(g.ii)){
  ## how large is the gap
  k  <- g.ii[j] - g.ii[j-1]
  if(k == 1){
    ## now set the labels
    text(p.ii[g.ii[j]], pm[4]+strheight("B")*.5, labels = tmp$top.gene[g.ii[j]], cex=.3, srt=90, xpd=NA,
         pos=4, offset = 0, font=3)
    ## add arrow
    arrows(p.ii[g.ii[j]], pm[4]+strheight("B")*.5, p.ii[g.ii[j]],
           pm[4], lwd=.3, length=0, xpd=NA)
    # arrows(p.ii[g.ii[j]], pm[4], p.ii[g.ii[j]],
    #        tmp$id[g.ii[j]]+2, lwd=.3, length=0, xpd=NA, lty=2)
  }else{
    ## define mean of the points in the gap
    m  <- mean(p.ii[(g.ii[j-1]+1):g.ii[j]])
    ## length of the whole block
    bw <- k*strwidth("B")*.5
    
    print(tmp$candidate.gene[(g.ii[j-1]+1):g.ii[j]])
    print(p.ii[(g.ii[j-1]+1):g.ii[j]])
    print(m)
    print(seq(m-(bw/2), m+(bw/2), length.out=k))
    
    ## now set the labels
    text(seq(m-(bw/2), m+(bw/2), length.out=k), pm[4]+strheight("B")*.5, 
         labels = tmp$top.gene[(g.ii[j-1]+1):g.ii[j]], cex=.3, srt=90, xpd=NA,
         pos=4, offset = 0, font=3)
    ## add arrows
    arrows(seq(m-(bw/2), m+(bw/2), length.out=k), pm[4]+strheight("B")*.5, p.ii[(g.ii[j-1]+1):g.ii[j]],
           pm[4], lwd=.3, length=0, xpd=NA)
    # arrows(p.ii[(g.ii[j-1]+1):g.ii[j]], pm[4], p.ii[(g.ii[j-1]+1):g.ii[j]],
    #        tmp$id[(g.ii[j-1]+1):g.ii[j]]+2, lwd=.3, length=0, xpd=NA, lty=2)
  }
}

## add categories on top
points(region.count$x.pos, ifelse(region.count$n.assoc >= 50, 50, region.count$n.assoc), pch=21, cex=region.count$n.categories/15,
       bg="white", lwd=.3, xpd=NA)


#----------------------------------#
##-- lower part for single vars --##
#----------------------------------#

## adapt plotting parameters
par(mar=c(3.5,1.5,.1,.5))

## empty plot
plot(c(chr.dat$start[1], chr.dat$end[23]), c(0,max(res.all$y.pos.v2)), xaxt="n", yaxt="n",
     xlab="Chromosomal position", ylab="", type="n", ylim=rev(c(0,max(res.all$y.pos.v2))))
## add axis
axis(1, lwd=.5, at=chr.dat$mid, labels=1:23)
## get plotting coordinates
pm <- par("usr")
## add phenotype categories
rect(pm[1], cl.dat2$start, pm[2], cl.dat2$end, col=colorspace::lighten(cl.dat2$cl, .6), border=cl.dat2$cl, lwd=.1)
abline(v=chr.dat$start, lwd=.3, lty=2, col="white")
## add points
tmp <- res.all[ unreported != "unreported"]
# points(tmp$x.pos, tmp$y.pos.v2, pch=21, bg=tmp$cl.rep, cex=.25, lwd=.1, xpd=NA)
points(tmp$x.pos, tmp$y.pos.v2, pch=21, bg="grey80", cex=.25, lwd=.1, xpd=NA)
tmp <- res.all[ unreported == "unreported"]
# points(tmp$x.pos, tmp$y.pos.v2, pch=21, bg=tmp$cl.rep, cex=.25, lwd=.1, xpd=NA)
points(tmp$x.pos, tmp$y.pos.v2, pch=21, bg="black", cex=.25, lwd=.1, xpd=NA)

## add line to devide
abline(h=721, lwd=.5)
## add axis labels
text(pm[1]-(pm[2]-pm[1])*.01, 721+(724/2), labels = "Disease elsewhere", cex=.5, xpd=NA, srt=90, pos=3, offset=0)
text(pm[1]-(pm[2]-pm[1])*.01, 721/2, labels = "Disease mostly seen\nin primary care", cex=.5, xpd=NA, srt=90, pos=3, offset=0)

## add legend
legend(pm[1]-(pm[2]-pm[1])*.02, pm[3]-(pm[4]-pm[3])*.1, lty=0, pch=22, pt.lwd=.1, bty="n", xpd=NA,
       pt.bg=c(cl.dat$cl, "black", "grey80"), legend = c(cl.dat$category, "unreported", "reported"),
       cex=.4, ncol=8, pt.cex=1)

dev.off()

########################################################################
########################################################################
####                   Revision NatGen - 26/09/2023                 ####
########################################################################
########################################################################


#########################################
####  candidate gene score vs OT-V2G ####
#########################################

## define number of times, for which v2g top genes that conflicts with
## top, robust coloc gene
tmp <- res.all

## create new column to ease matching
tmp[, v2g.gene := sapply(v2g.top, function(x) strsplit(x, "_")[[1]][1])]
## same for coloc (only use strong justification, PP.H4 > .8)
tmp[, coloc.gene := apply(tmp[, c("eQTL_summary.coloc", "pQTL_summary.coloc"), with=F], 1, function(x){
  
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
    # e.gen$PP     <- ifelse(e.gen$PP > .8, 1, ifelse(e.gen$PP > .5, .5, 0))
  }else{
    e.gen        <- data.frame(gene="", PP=NA)
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
    # e.prot$PP     <- ifelse(e.prot$PP > .8, 1, ifelse(e.prot$PP > .5, .5, 0))
  }else{
    e.prot        <- data.frame(gene="", PP=NA)
  }
  
  ## combine and report highest
  tmp.gen <- na.omit(unique(rbind(e.gen, e.prot)))
  print(tmp.gen)
  
  ## report if any
  if(nrow(subset(tmp.gen, PP >= .8)) > 0){
    ## order 
    tmp.gen <- tmp.gen[ order(-tmp.gen$PP),]
    return(tmp.gen$gene[1])
  }else{
    return(NA)
  }

  
})]

## how many with different annotations
nrow(tmp[ !is.na(v2g.gene) & !is.na(coloc.gene) & v2g.gene != coloc.gene])

############################################
####            heritability            ####
############################################

## look at heritability vs samples size
cor.test(res.ldsc$h2_est, res.ldsc$cases, m="s", u="p")$p.value

## compare heritability across groups
boxplot(h2_est ~ category, res.ldsc)

## simple test
summary(lm(h2_est ~ cases + category, res.ldsc))

#                                     Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)                      3.936e-03  5.299e-04   7.428 1.89e-13 ***
#   cases                            3.592e-07  9.124e-09  39.373  < 2e-16 ***
#   categorycongenital anomalies    -2.964e-03  1.023e-03  -2.897 0.003830 ** 
#   categorydermatologic            -3.012e-03  8.383e-04  -3.593 0.000338 ***
#   categorydigestive               -1.789e-03  7.366e-04  -2.428 0.015290 *  
#   categoryendocrine/metabolic     -5.871e-04  7.755e-04  -0.757 0.449107    
#   categorygenitourinary           -1.143e-03  7.288e-04  -1.568 0.117182    
#   categoryhematopoietic           -3.534e-03  1.096e-03  -3.223 0.001298 ** 
#   categoryinfectious diseases     -4.694e-03  9.778e-04  -4.801 1.75e-06 ***
#   categoryinjuries & poisonings   -3.344e-03  8.308e-04  -4.025 5.99e-05 ***
#   categorymental disorders        -1.712e-03  9.160e-04  -1.869 0.061887 .  
#   categorymusculoskeletal         -1.962e-03  8.008e-04  -2.450 0.014422 *  
#   categoryneoplasms               -1.856e-03  7.612e-04  -2.438 0.014893 *  
#   categoryneurological            -2.430e-03  8.854e-04  -2.745 0.006133 ** 
#   categoryother                   -6.573e-03  1.594e-03  -4.123 3.95e-05 ***
#   categorypregnancy complications -2.796e-03  1.086e-03  -2.574 0.010165 *  
#   categoryrespiratory             -3.255e-03  8.930e-04  -3.645 0.000277 ***
#   categorysense organs            -2.483e-03  7.848e-04  -3.164 0.001589 ** 
#   categorysymptoms                -5.842e-03  1.212e-03  -4.819 1.60e-06 ***






