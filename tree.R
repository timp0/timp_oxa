library(tidyverse)
library(googlesheets)

workdir="/atium/Data/Nanopore/Analysis/170104_reassemble/annot/"
outdir="~/Dropbox/Data/Nanopore/170104_reassemble"

refdir="/atium/Data/Reference/kpneumo"

gs_auth(token = "~/googlesheets_token.rds")

make.tree = function(workdir, outdir, label="nanocontigs", assembly="nanopore.canu", kpneumo.refs, parsnp.label="nano_parsnp", full.tree=F) {

    ##Raw nanopore assembly tree
    ##copy nanopore assemblies
    contig.dir=file.path(workdir, label)
    if (!dir.exists(contig.dir)) {
        dir.create(contig.dir)
    }
    contig=dataloc %>%
        rename_(afile=assembly) %>%
        filter(!is.na(afile)) %>%
        mutate(dest.file=file.path(contig.dir, paste0("isolate_", trish.id, ".fasta"))) %>%
        select(afile, dest.file)
    file.copy(contig$afile, contig$dest.file)       
    
    ##copy reference genome
    if (full.tree) {
        kpneumo.refs=kpneumo.refs %>%
            filter(Strain!="NA") %>%
            mutate(Strain=gsub(" ", "", Strain)) %>%
            mutate(Strain=gsub("[[:punct:]]", "", Strain))
        file.copy(kpneumo.refs$dest.name, file.path(contig.dir, paste0(kpneumo.refs$Strain, ".fasta.gz")))
        system(paste0("gunzip -f ", file.path(contig.dir, "*.gz")))
    } else {
        ##Just HS11286 for now
        file.copy(kpneumo.refs$dest.name[kpneumo.refs$Strain=="HS11286"&(!is.na(kpneumo.refs$Strain))], file.path(contig.dir, "HS11286.fasta.gz"))
        ##Got to ungzip these for parsnp (obnoxiously)
        system(paste0("gunzip -f ", file.path(contig.dir, "HS11286.fasta.gz")))
    }
    ##parsnp
    system(paste0("~/Code/Harvest-Linux64-v1.1.2/parsnp -r ", file.path(contig.dir, "HS11286.fasta"), " -d ", contig.dir, " -p 10 -o ./", parsnp.label))
    
    ##Count snps
    vcfloc=file.path(workdir,parsnp.label, "parsnp.vcf")  
    system(paste0("~/Code/Harvest-Linux64-v1.1.2/harvesttools -i ", file.path(workdir, parsnp.label, "parsnp.ggr"), " -V ", vcfloc))
    count.snps(vcfloc)
    
    ##copy tree file out 
    system(paste0("cp -R ", file.path(workdir, parsnp.label), " ", file.path(outdir)))
}


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

#ref.cols=c("assembly_accession", "bioproject", "biosample", "wgs_master", "refseq_category", "taxid", "species_taxid", "organism_name",
#           "infraspecific_name", "isolate", "version_status", "assembly_level", "release_type", "genome_rep", "seq_rel_date", "asm_name",
#           "submitter", "gbrs_paired_asm", "paired_asm_comp", "ftp_path", "excluded_from_refseq")

#kpneumo.refs=read_tsv(file.path(outdir, "assembly_summary.txt"), comment="#")

##All kpneumo refs
kpneumo.refs=read_csv(file.path(outdir, "genomes_proks.csv")) %>%
    rename(org.name=`#Organism/Name`) %>%
    mutate(rootname=sapply(strsplit(`GenBank FTP`, split="/"), function(x) {x[length(x)]})) %>%
    mutate(source.name=paste0(`GenBank FTP`, "/", rootname, "_genomic.fna.gz")) %>%
    mutate(dest.name=file.path(refdir, paste0(rootname, "_genomic.fna.gz")))

##Download 
if (FALSE) {

    for (i in 1:dim(kpneumo.refs)[1]) {
        ##Get all refs
        download.file(kpneumo.refs$source.name[i], destfile=kpneumo.refs$dest.name[i])
    }
    
}

setwd(workdir)

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1_WT3RQSVGvR97-asIHtIy0WWg49rAFiWWQe9g8BWNiQ/edit?usp=sharing")
dataloc=gs_read(fullsheet, ws="KPneumo0123")


##Raw nanopore assembly tree

make.tree(workdir, outdir, label="nanocontigs", assembly="nanopore.canu", kpneumo.refs, parsnp.label="nano_parsnp")
make.tree(workdir, outdir, label="piloncontigs", assembly="pilon", kpneumo.refs, parsnp.label="pilon_parsnp")
make.tree(workdir, outdir, label="spadescontigs", assembly="illumina.spades", kpneumo.refs, parsnp.label="spades_parsnp")

##Next
##Generate massive tree with nanopore assemblies and all K. Pnemuo refs

make.tree(workdir, outdir, label="pilonfull", assembly="pilon", kpneumo.refs, parsnp.label="pilon_full_parsnp", full.tree=T)

