#!/bin/bash

#command I used for prokka (example with 10A using Illumina)

/home/shao4/prokka/bin/prokka --outdir /atium/Data/NGS/Aligned/161027_SSSSS/SPAdes/assemblies/10A_genus --genus Klebsiella --usegenus --addgenes --prefix 10A_genus /atium/Data/NGS/Aligned/161027_SSSSS/SPAdes/spade_contigs/Illumina/10A_contigs.fasta

#note that on some of them I added the --addgenes flag - I didn't find it very helpful, as it seemed to just duplicate your results

#with nanopore
/home/shao4/prokka/bin/prokka --outdir /atium/Data/NGS/Aligned/161027_SSSSS/SPAdes/assemblies/1A_nanopore_genus --genus Klebsiella --usegenus --prefix 1A_nanopore /atium/Data/NGS/Aligned/161027_SSSSS/SPAdes/spade_contigs/nanopore_Illumina/Nanopore/161109_SSSSS_1A.nanopore.contigs.fasta

#building the custom database for Klebsiella (notes also in evernote):
# first, download all of the genbank files from NCBI for various K. pneumo strains - I have them downloaded at /atium/Data/NGS/Aligned/161027_SSSSS/SSSSS_ref/
#follow instructions as listed on website

/home/shao4/prokka/bin/prokka-genbank_to_fasta_db *.gbk > Klebsiella.faa
cd-hit -i Klebsiella.faa -o Klebsiella -T 0 -M 0 -g 1 -s 0.8 -c 0.9
rm -fv Klebsiella.clstr Klebsiella.faa
makeblastdb -dbtype prot -in Klebsiella
mv Klebsiella.p* /home/shao4/prokka/db/genus



