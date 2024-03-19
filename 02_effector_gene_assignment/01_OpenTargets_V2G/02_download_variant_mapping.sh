#!/bin/sh

## script to download variant mapping files from OpenTargets
## Maik Pietzner                                  07/01/2022

#BSUB -e <path to file>
#BSUB -o <path to file>

## force LSF to allocate all cores on one node; implies
## that no job can be submitted if it tries to allocate more cores than a single
## node can offer
#BSUB -R "span[hosts=1]" 

## Request 2GB of memory/core
#BSUB -R "rusage[mem=2GB]"

## how many threads
#BSUB -n 2

## name of the job (do as an array per chromosome)
#BSUB -J op_download[199,136,150,145,175,118,151]%5

## to which queue to be send (type bqueues to see all)
#BSUB -q normal

## get the output to be collated
echo "Job ID: $LSB_JOBINDEX"
# export num=$((${LSB_JOBINDEX} - 1))
export num=${LSB_JOBINDEX}

## change to download destination
cd <path to file>

## create query depending on numeric identifier
if [ ${num} -lt 10 ]; then

## delete possible old file
rm part-0000${num}-b7926934-3aa8-4f76-a97f-488a1eaa7262-c000.json
wget http://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/20022712/lut/variant-index/part-0000${num}-b7926934-3aa8-4f76-a97f-488a1eaa7262-c000.json

elif [[ ${num} -ge 10 && ${num} -lt 100 ]]; then

## delete possible old file
rm part-000${num}-b7926934-3aa8-4f76-a97f-488a1eaa7262-c000.json
wget http://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/20022712/lut/variant-index/part-000${num}-b7926934-3aa8-4f76-a97f-488a1eaa7262-c000.json

else

## delete possible old file
rm part-00${num}-b7926934-3aa8-4f76-a97f-488a1eaa7262-c000.json
wget http://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/20022712/lut/variant-index/part-00${num}-b7926934-3aa8-4f76-a97f-488a1eaa7262-c000.json

fi