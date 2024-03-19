#!/usr/bin/env Rscript

## script to run cis-eQTL colocalisation at GWAS regions
## Maik Pietzner
rm(list=ls())

## get the arguments from the command line
args <- commandArgs(trailingOnly=T)

## little options
options(stringsAsFactors = F)
## avoid conversion of numbers
options(scipen = 1)
# print(R.Version())

setwd("<path to file>")

## --> import parameters <-- ##

## get regional coordinates
pheno <- args[1]
chr.s <- as.numeric(args[2])
pos.l <- as.numeric(args[3])
pos.s <- as.numeric(args[4])
pos.e <- as.numeric(args[5])

## get regional coordinates
# tmp   <- read.table("input/input.eQTL.pipeline.txt")
# j     <- 5778
# pheno <- tmp$V1[j]
# chr.s <- tmp$V2[j]
# pos.l <- tmp$V3[j]
# pos.s <- tmp$V4[j]
# pos.e <- tmp$V5[j]

#-----------------------------------------#
##--     load regional assoc stats     --##
#-----------------------------------------#

## package for faster data handling
require(data.table)

## read the relevant data
res            <- paste0("zcat /sc-projects/sc-proj-computational-medicine/people/Maik/14_phecode_GWAS/01_GWAS/gwas_results/", pheno, ".allchr.results.gz | awk -v chr=", chr.s, " -v low=", pos.s, " -v upp=", pos.e, 
                         " '{if(($1 == chr && $2 >= low && $2 <= upp) || NR == 1) print $0}'")
res            <- data.table::fread(cmd = res, sep = " ", header = T, data.table = F)

## drop SNPs that have possibly failed in REGENIE
res            <- subset(res, is.na(EXTRA))

## create MarkerName to enable mapping to the LD file
res$MarkerName <- apply(res, 1, function(x){
  paste0("chr", as.numeric(x[1]), ":", as.numeric(x[2]), "_", paste(sort(x[4:5]), collapse = "_"))
})

#-----------------------------------------#
##--        lift over to GRCh38        --##
#-----------------------------------------#

## packages needed
require(biomaRt)
require(liftOver)

# ## change chromosome if needed
# if(chr.s == 23){
#   chr.s <- "X"
# }

## import chain
chain      <- import.chain("input/hg19ToHg38.over.chain")

## create GRanges object
grObject   <- GRanges(seqnames = paste0("chr", ifelse(chr.s == 23, "X", chr.s)), ranges=IRanges(start = res$GENPOS, end = res$GENPOS, names=res$MarkerName))

## now do the mapping
tmp        <- as.data.frame(liftOver(grObject, chain))
## rename
names(tmp) <- c("group", "MarkerName", "seqnames", "pos.hg38", "end", "width", "strand")
## add to the data
res        <- merge(res, tmp[, c("MarkerName", "pos.hg38")])

#-----------------------------------------#
##--       identify genes closeby      --##
#-----------------------------------------#

## get data on build 37
gene.ensembl   <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37) 

## get all genes in the region (adapt a larger region, just to be sure to get any output)
tmp.genes      <- getBM(attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype'),
                        filters = c('chromosome_name','start','end'),
                        values = list(ifelse(chr.s == 23, "X", chr.s), pos.s-5e5, pos.e+5e5),
                        mart = gene.ensembl)
## restrict to protein encoding genes for now
tmp.genes      <- subset(tmp.genes, gene_biotype %in% c("protein_coding", "processed_transcript"))

## subset to genes no more than 500kb away from the lead variant (gene body not TSS)
tmp.genes$dist <- apply(tmp.genes[, c("start_position", "end_position")], 1, function(x) min(abs(x-pos.l)))
tmp.genes      <- subset(tmp.genes, dist <= 5e5)

#-----------------------------------------#
##--          start coloc              --##
#-----------------------------------------#

## do only if at least one gene has been found
if(nrow(tmp.genes) > 0){
  
  #-----------------------------------------#
  ##--          import LD matrix         --##
  #-----------------------------------------#
  
  ## write list of SNPs to be queried to file 
  write.table(res$ID, paste("tmpdir/snplist", pheno, chr.s, pos.s, pos.e, "lst", sep="."), col.names = F, row.names = F, quote = F)
  
  ## snp data to be queried
  tmp.z        <- res[, c("ID", "CHROM", "GENPOS", "ALLELE0", "ALLELE1")]
  names(tmp.z) <- c("rsid", "chromosome", "position", "allele1", "allele2")
  
  ## adopt chromosome if needed
  if(chr.s < 10){
    print("tada")
    tmp.z$chromosome <- paste0("0", tmp.z$chromosome)
  }else if(chr.s == 23){
    tmp.z$chromosome <- "X"
  }
  
  ## check for input
  print(head(tmp.z))
  
  ## write to file
  write.table(tmp.z, paste("tmpdir/snpz", pheno, chr.s, pos.s, pos.e, "z", sep="."), row.names = F, quote = F)
  
  ## --> create master file for LDstore2 <-- ##
  
  ## assign entries
  m.file      <- data.frame(z=paste("tmpdir/snpz", pheno, chr.s, pos.s, pos.e, "z", sep="."),
                            bgen=paste("tmpdir/filtered", pheno, chr.s, pos.s, pos.e, "bgen", sep="."),
                            bgi=paste("tmpdir/filtered", pheno, chr.s, pos.s, pos.e, "bgen.bgi", sep="."),
                            ld=paste("tmpdir/ld", pheno, chr.s, pos.s, pos.e, "ld", sep="."),
                            incl="<path to file>",
                            n_samples=ifelse(chr.s != 23, 487409, 486757))
  
  ## write to file
  write.table(m.file, paste("tmpdir/master", pheno, chr.s, pos.s, pos.e, "z", sep="."), sep=";", row.names = F, quote = F)
  
  cat("--------------------------------------\n")
  cat("computing LD matrix\n")
  
  ## submit the job
  system(paste("./scripts/get_LD_matrix_ldstore.sh", chr.s, pos.s, pos.e, pheno))
  
  ## read in matrix
  ld         <- fread(paste("tmpdir/ld", pheno, chr.s, pos.s, pos.e, "ld", sep="."), data.table = F)
  ## print for checking
  print(ld[1:5,1:5])
  
  ## create identifier column in results to keep the mapping to the LD matrix
  res$snp.id <- 1:nrow(res) 
  
  cat("Done\n")
  cat("--------------------------------------\n")
  
  #-----------------------------------------#
  ##--     start coloc for each gene     --##
  #-----------------------------------------#
  
  ## import function to do so
  source("scripts/coloc_gene_gtex_tissues.R")
  
  ## run across all genes
  res.genes <- lapply(tmp.genes$ensembl_gene_id, function(g){
    
    print(g)
    
    ## run coloc for specific gene
    res.tmp                 <- coloc.eQTL(res, chr.s, pos.s, pos.e, g, ld)
    
    print(head(res.tmp))
    
    ## return only if any
    if(!is.null(res.tmp)){
      ## add gene
      res.tmp$ensembl_gene_id <- g
      ## add readable name
      res.tmp$gene_name       <- tmp.genes$external_gene_name[which(tmp.genes$ensembl_gene_id == g)]
      
      ## return results
      return(res.tmp)
    }
  })
  
  ## collate and combine
  res.genes <- do.call(rbind, res.genes)
  
  ## write results to file
  write.table(res.genes, paste("output/eQTL", pheno, chr.s, pos.s, pos.e, "txt", sep = "."), sep="\t", row.names=F)
  
}else{
  
  cat("no protein-encoding genes found closeby\n")
  ## export empty results file to keep track on whether the job was run successfully
  write.table(data.frame(tissue=NA), paste("output/eQTL", pheno, chr.s, pos.s, pos.e, "txt", sep = "."), sep="\t", row.names=F)
  
}
