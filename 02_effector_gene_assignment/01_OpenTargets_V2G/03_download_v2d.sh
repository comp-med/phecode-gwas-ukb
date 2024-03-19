#!/bin/sh

## script to download v2d files from OpenTargets
## Maik Pietzner                                  12/01/2022

#BSUB -e <path to file>
#BSUB -o <path to file>

## force LSF to allocate all cores on one node; implies
## that no job can be submitted if it tries to allocate more cores than a single
## node can offer
#BSUB -R "span[hosts=1]" 

## Request 2GB of memory/core
#BSUB -R "rusage[mem=2GB]"

## how many threads
#BSUB -n 4

## name of the job (do as an array per chromosome)
#BSUB -J op_download[3-5]%3

## to which queue to be send (type bqueues to see all)
#BSUB -q normal

## get the output to be collated
echo "Job ID: $LSB_JOBINDEX"
export num=$((${LSB_JOBINDEX} - 1))
# export num=${LSB_JOBINDEX}

## change to download destination
cd <path to file>

echo ${num}

## create query depending on numeric identifier
if [ ${num} -lt 10 ]; then

## delete possible old file
rm -rf part-0000${num}-bd4a8606-3a8a-4862-bd1d-1828cc2d6332-c000.json
wget http://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/latest/v2d/part-0000${num}-bd4a8606-3a8a-4862-bd1d-1828cc2d6332-c000.json

elif [[ ${num} -ge 10 && ${num} -lt 100 ]]; then

## delete possible old file
rm -rf part-000${num}-bd4a8606-3a8a-4862-bd1d-1828cc2d6332-c000.json
wget http://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/latest/v2d/part-000${num}-bd4a8606-3a8a-4862-bd1d-1828cc2d6332-c000.json

else

## delete possible old file
rm -rf part-00${num}-bd4a8606-3a8a-4862-bd1d-1828cc2d6332-c000.json
wget http://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/latest/v2d/part-00${num}-bd4a8606-3a8a-4862-bd1d-1828cc2d6332-c000.json

fi