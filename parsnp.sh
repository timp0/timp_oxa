#parsnp scripts

#for Illumina
/home/shao4/Harvest-Linux64-v1.1.2/parsnp -c -r /atium/Data/NGS/Aligned/161027_SSSSS/SPAdes/spade_contigs/nanopore_Illumina/Illumina/HS11286.fasta -d /atium/Data/NGS/Aligned/161027_SSSSS/SPAdes/spade_contigs/nanopore_Illumina/Illumina/ -p 10 -o ./Illumina_parsnp

#for Nanopore
/home/shao4/Harvest-Linux64-v1.1.2/parsnp -c -r /atium/Data/NGS/Aligned/161027_SSSSS/SPAdes/spade_contigs/nanopore_Illumina/Nanopore/HS11286.fasta -d /atium/Data/NGS/Aligned/161027_SSSSS/SPAdes/spade_contigs/nanopore_Illumina/Nanopore/ -p 10 -o ./Nanopore_parsnp


#the -c tag was left on because the first time I ran it a sample or two didn't get put on the tree for whatever reason (-c stands for curated genome directory); -r is your reference file, -d is the directory where everything you are making into a tree is, -p is number of threads, -o is output directory
