#!/bin/bash

#makeblastdb -in nucleotide_fasta_protein_homolog_model.fasta -out CARDblast -dbtype nucl -input_type fasta
#lastdb -uNEAR -R01 CARDlast nucleotide_fasta_protein_homolog_model.fasta

assembly=811A.nanopore.contigs.fasta
##takes just before .
label="${assembly%%.*}"


##Blast reads against CARD
##Old per read code
#cat 160223_VRE7_recall_2dhq.fa | parallel --block 100k --recstart '>' --pipe \
#    blastn -db CARDblast -query - -gapopen 1 -gapextend 2 -word_size 9 -reward 1 -evalue 10 \
#    -outfmt '"6 qseqid sseqid pident length qlen slen evalue bitscore"' > nanoblast2

##Keep just bet hit for region - a lot of CARD is very similar, causing problems if I don't set small best_hit_score_edge
##See https://www.ncbi.nlm.nih.gov/books/NBK279668/#usermanual.BestHits_filtering_algorithm for info
blastn -db /atium/Data/Reference/CARD/CARDblast -query ${assembly} -evalue 1e-50 -outfmt 6 -best_hit_overhang .1 -best_hit_score_edge 1e-5 >${label}.cardblast.tsv


##RGI
python ~/Code/rgi_card/rgi.py -t contig -i ${assembly} -n 8 -o ${label}.rgicard.json
python ~/Code/rgi_card/convertJsonToTSV.py -i ${label}.rgicard.json -o ${label}.rgicard


##makeblastdb -in isecp1.fasta -out isecp1blast -dbtype nucl -input_type fasta
##makeblastdb -in mutation_dbase.fasta -out mutdbase -dbtype nucl -input_type fasta
blastn -db genedb/mutdbase -query ${assembly} -outfmt 6 >${label}.mutblast.tsv
blastn -task blastn-short -db genedb/isecp1blast -query ${assembly} -outfmt 6 >${label}.iseblast.tsv



##From other code
#last-train -P10 CARDlast 160223_VRE7_recall_2dhq.fa  > CARD.par
#cp /atium/Data/Nanopore/Analysis/161202_last/r73.par .
#lastal -P10 -p r73.par CARDlast 160223_VRE7_recall_2dhq.fa | last-split -m1e-6 > lastcard.maf

