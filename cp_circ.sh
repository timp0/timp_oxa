#!/bin/bash


for i in 1 2 4 6 7 8 9 10 12 ;
do
    tiggfa=/atium/Data/Nanopore/Analysis/170104_reassemble/canu/$i.nanopore/$i.nanopore.contigs.gfa
    cp $tiggfa ~/Dropbox/Lab/oxa/canu_gfa
done
