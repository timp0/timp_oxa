#!/bin/bash

##bash blast_reads.sh fqpath generefpath

datadir=/atium/Data/Nanopore/Analysis/170104_reassemble/blast/samp7
outdir=~/Dropbox/Lab/oxa

for i in $datadir/*fastq.gz ;
do
    samp=`basename $i .fastq.gz`
    readdir=`dirname $i`

    ##make fa for blast if not already done
    if [ ! -f $readdir/$samp.fasta ] ; then
	seqtk seq -a $i > $readdir/$samp.fasta
    fi
done

for i in $datadir/*.fasta ;
do

    samp=`basename $i .fasta`
    readdir=`dirname $i`
    
    ##create blast db in the dir of the reads if not already there
    if ls $readdir/$samp.db* 1> /dev/null 2>&1; then
	echo 'db already exists'
    else
	makeblastdb -in $readdir/$samp.fasta -out $readdir/$samp.db -dbtype nucl
    fi
    
    
    for j in /atium/Data/Nanopore/Analysis/170104_reassemble/oxa_refs/*fasta ;
    do
	##search
	gene=`basename $j .fasta`
	##only need tsv for this part, I think
	blastn -query $j -db $readdir/$samp.db -outfmt 7 -out $outdir/$samp.$gene.tsv
    done
    
done
    


