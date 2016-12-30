#!/bin/bash


##From Steph (altered by WT)
##Making K.pneumo dbase
~/Code/prokka/bin/prokka-genbank_to_fasta_db /atium/Data/NGS/Aligned/161027_SSSSS/SSSSS_ref/*.gbk > Klebsiella.faa

#building the custom database for Klebsiella (notes also in evernote):
# first, download all of the genbank files from NCBI for various K. pneumo strains - I have them downloaded at /atium/Data/NGS/Aligned/161027_SSSSS/SSSSS_ref/
#follow instructions as listed on website

cd-hit -i Klebsiella.faa -o Klebsiella -T 8 -M 24000 -g 1 -s 0.8 -c 0.9
rm -fv Klebsiella.clstr Klebsiella.faa
makeblastdb -dbtype prot -in Klebsiella
mv Klebsiella.p* ~/Code/prokka/db/genus


##From Steph (altered by WT)
#command I used for prokka 
label=811A

~/Code/prokka/bin/prokka --outdir /atium/Data/Nanopore/Analysis/161223_kpneumo/${label}_genus --genus Klebsiella --usegenus \
    --prefix 811A_genus /atium/Data/NGS/Aligned/161027_SSSSS/SPAdes/spade_contigs/nanopore_Illumina/Nanopore/811A.nanopore.contigs.fasta

#note that on some of them I added the --addgenes flag - I didn't find it very helpful, as it seemed to just duplicate your results




