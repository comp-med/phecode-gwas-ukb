#!/bin/sh

## script to create a LD-matrix
## Maik Pietzner 02/03/2022

## export location of files
export dir=<path to file>/bgen_files

## get the chromosome
export chr=${1}
export lowpos=${2}
export uppos=${3}
export pheno=${4}

## cd in the correct directory
cd <path to file>

echo "Chromosome ${chr} : Locus start ${lowpos} : Locus end ${uppos}"

if [ $chr -eq 23 ]; then

  ## create subset BGEN file (rsids of interest)
  <path to file>/bgenix \
  -g <path to file>/bgen_files/ukb_imp_chrX_v3.bgen \
  -incl-rsids tmpdir/snplist.${pheno}.${chr}.${lowpos}.${uppos}.lst > tmpdir/filtered.${pheno}.${chr}.${lowpos}.${uppos}.bgen

else

  ## create subset BGEN file (rsids of interest)
  <path to file>/bgenix \
  -g <path to file>/bgen_files/ukb_imp_chr${chr}_v3.bgen \
  -incl-rsids tmpdir/snplist.${pheno}.${chr}.${lowpos}.${uppos}.lst > tmpdir/filtered.${pheno}.${chr}.${lowpos}.${uppos}.bgen
  
fi

## create corresponding index file
<path to file>/bgenix \
-g tmpdir/filtered.${pheno}.${chr}.${lowpos}.${uppos}.bgen \
-index


## run LD correlation
<path to file>/ldstore_v2.0_x86_64/ldstore_v2.0_x86_64 \
--in-files tmpdir/master.${pheno}.${chr}.${lowpos}.${uppos}.z \
--write-text \
--n-threads 30 \
--read-only-bgen

## remove BGEN no longer needed
rm tmpdir/filtered.${pheno}.${chr}.${lowpos}.${uppos}.bgen
rm tmpdir/filtered.${pheno}.${chr}.${lowpos}.${uppos}.bgen.bgi
rm tmpdir/master.${pheno}.${chr}.${lowpos}.${uppos}.z
rm tmpdir/snpz.${pheno}.${chr}.${lowpos}.${uppos}.z