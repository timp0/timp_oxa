library(tidyverse)

##Want to use gsheets directly, but right now doesn't work easily, skip for now
#library(googlesheets)
#options(httr_oob_default=T)
#gs_auth(new_user=T)

count.snps = function(vcfloc) {
    ##Count snps in vcf file

    vcfs=read_tsv(vcfloc, comment="##")

    ##ok - for each sample column, need to look at delta for all other sample columns, then sum to get that square on the table

    ##Interesting names are names(vcfs)[10:length(vcfs)]

    comp=as_data_frame(t(combn(names(vcfs[10:length(vcfs)]),2))) %>%
        mutate(val=colSums(abs(vcfs[V1]-vcfs[V2])))

    comp=rbind(comp, data_frame(V1=comp$V2, V2=comp$V1, val=comp$val))

    write_csv(spread(comp, V2, val), paste0(vcfloc, ".csv"))
}



workdir="/atium/Data/Nanopore/Analysis/161223_kpneumo/trees"
outdir="~/Dropbox/Data/Nanopore/161230_oxatree"

refdir="/atium/Data/Reference/kpneumo"

#ref.cols=c("assembly_accession", "bioproject", "biosample", "wgs_master", "refseq_category", "taxid", "species_taxid", "organism_name",
#           "infraspecific_name", "isolate", "version_status", "assembly_level", "release_type", "genome_rep", "seq_rel_date", "asm_name",
#           "submitter", "gbrs_paired_asm", "paired_asm_comp", "ftp_path", "excluded_from_refseq")

#kpneumo.refs=read_tsv(file.path(outdir, "assembly_summary.txt"), comment="#")

kpneumo.refs=read_csv(file.path(outdir, "genomes_proks.csv")) %>%
    rename(org.name=`#Organism/Name`) %>%
    mutate(rootname=sapply(strsplit(`GenBank FTP`, split="/"), function(x) {x[length(x)]})) %>%
    mutate(source.name=paste0(`GenBank FTP`, "/", rootname, "_genomic.fna.gz")) %>%
    mutate(dest.name=file.path(refdir, paste0(rootname, "_genomic.fna.gz")))

if (FALSE) {

    for (i in 1:dim(kpneumo.refs)[1]) {
        ##Get all refs
        download.file(kpneumo.refs$source.name[i], destfile=kpneumo.refs$dest.name[i])
    }
    
}

setwd(workdir)

##Load CSV with sample info

dataloc=read_csv(file.path(outdir, "dataloc.csv"))

##copy nanopore assemblies

contig.dir=file.path(workdir, "nanocontigs")

if (!dir.exists(contig.dir)) {
    dir.create(contig.dir)
}

nano.contig=dataloc %>%
    filter(Sequenced=="both") %>%
    mutate(dest.file=file.path(contig.dir, paste0(trish.id, ".fasta"))) %>%
    select(canu.assembly, dest.file)

file.copy(nano.contig$canu.assembly, nano.contig$dest.file)       

##copy reference genome

##Just HS11286 for now

file.copy(kpneumo.refs$dest.name[kpneumo.refs$Strain=="HS11286"&(!is.na(kpneumo.refs$Strain))], file.path(contig.dir, "HS11286.fasta.gz"))
##Got to ungzip these for parsnp (obnoxiously)
system(paste0("gunzip ", file.path(contig.dir, "HS11286.fasta.gz")))

##parsnp

system(paste0("~/Code/Harvest-Linux64-v1.1.2/parsnp -r ", file.path(contig.dir, "HS11286.fasta"), " -d ", contig.dir, " -p 10 -o ./Nanopore_parsnp"))

##Count snps

vcfloc=file.path(workdir, "Nanopore_parsnp", "parsnp.vcf")
    
system(paste0("~/Code/Harvest-Linux64-v1.1.2/harvesttools -i ", file.path(workdir, "Nanopore_parsnp", "parsnp.ggr"), " -V ", vcfloc))

count.snps(vcfloc)

##copy tree file out 

system(paste0("cp -R ", file.path(workdir, "Nanopore_parsnp"), " ", file.path(outdir, "Nanopore_parsnp")))

##copy illumina assemblies

contig.dir=file.path(workdir, "illcontigs")

if (!dir.exists(contig.dir)) {
    dir.create(contig.dir)
}

ill.contig=dataloc %>%
    mutate(dest.file=file.path(contig.dir, paste0(trish.id, ".fasta"))) %>%
    select(spades.assembly, dest.file)

file.copy(ill.contig$spades.assembly, ill.contig$dest.file)       


##copy reference genome

file.copy(kpneumo.refs$dest.name[kpneumo.refs$Strain=="HS11286"&(!is.na(kpneumo.refs$Strain))], file.path(contig.dir, "HS11286.fasta.gz"))
##Got to ungzip these for parsnp (obnoxiously)
system(paste0("gunzip ", file.path(contig.dir, "HS11286.fasta.gz")))

##parsnp

system(paste0("~/Code/Harvest-Linux64-v1.1.2/parsnp -r ", file.path(contig.dir, "HS11286.fasta"), " -d ", contig.dir, " -p 10 -o ./Illumina_parsnp"))

##Count snps

vcfloc=file.path(workdir, "Illumina_parsnp", "parsnp.vcf")
    
system(paste0("~/Code/Harvest-Linux64-v1.1.2/harvesttools -i ", file.path(workdir, "Illumina_parsnp", "parsnp.ggr"), " -V ", vcfloc))

count.snps(vcfloc)

##copy tree file out 

system(paste0("cp -R ", file.path(workdir, "Illumina_parsnp"), " ", file.path(outdir, "Illumina_parsnp")))

##Next
##Generate massive tree with nanopore assemblies and all K. Pnemuo refs
