#!/bin/bash

for namey in 161109_SSSSS_12A 
 
do

    filey=`ls ${namey}_*d.tar_2dhq.fastq.gz`
    echo $filey


    /home/shao4/canu-1.3/Linux-amd64/bin/canu -p ${namey}.nanopore -d assembly/${namey}.nanopore -genomeSize=5m -nanopore-raw ${filey}

done
