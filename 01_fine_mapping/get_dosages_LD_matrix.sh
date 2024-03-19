#!/bin/sh

## script to obtain SNP dosages to create an LD-matrix
## Maik Pietzner 19/11/2021

## export location of files
export dir=<path to file>/bgen_files

## get the chromosome
export chr=${1}
export lowpos=${2}
export uppos=${3}
export pheno=${4}

echo "Chromosome ${chr} : Locus start ${lowpos} : Locus end ${uppos}"

if [ ${chr} -eq 23 ]; then

  ## create subset bgen file get only SNPs in the data set
  <path to file>/bgenix \
  -g ${dir}/ukb_imp_chrX_v3.bgen \
  -incl-rsids tmpdir/snplist.${pheno}.${chr}.${lowpos}.${uppos}.lst > tmpdir/${pheno}.${chr}.${lowpos}.${uppos}.bgen
  
  ## create dosage file 
  <path to file>/qctool_v2.0.7 \
  -g tmpdir/${pheno}.${chr}.${lowpos}.${uppos}.bgen \
  -s ${dir}/ukb_imp_chrX_v3.sample \
  -incl-samples <path to file>.sample \
  -og - \
  -ofiletype dosage > tmpdir/tmp.${pheno}.${chr}.${lowpos}.${uppos}.dosage
  
else

  ## create subset bgen file get only SNPs in the data set
  <path to file>/bgenix \
  -g ${dir}/ukb_imp_chr${chr}_v3.bgen \
  -incl-rsids tmpdir/snplist.${pheno}.${chr}.${lowpos}.${uppos}.lst > tmpdir/${pheno}.${chr}.${lowpos}.${uppos}.bgen
  
  ## create dosage file 
  /sc-projects/sc-proj-computational-medicine/programs/qctool/build/release/qctool_v2.0.7 \
  -g tmpdir/${pheno}.${chr}.${lowpos}.${uppos}.bgen \
  -s ${dir}/ukb_imp_chr${chr}_v3.sample \
  -incl-samples <path to file>.sample \
  -og - \
  -ofiletype dosage > tmpdir/tmp.${pheno}.${chr}.${lowpos}.${uppos}.dosage

fi

rm tmpdir/${pheno}.${chr}.${lowpos}.${uppos}.bgen
rm tmpdir/snplist.${pheno}.${chr}.${lowpos}.${uppos}.lst
