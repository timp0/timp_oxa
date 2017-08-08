library(tidyverse)
library(googlesheets)

workdir="/atium/Data/Nanopore/Analysis/170516_oxapilon/annot/"
outdir="~/Dropbox/Data/Nanopore/170718_oxa"

refdir="/atium/Data/Reference/kpneumo"

gs_auth(token = "~/googlesheets_token.rds")


make.tree = function(dataloc, workdir, outdir, label="nanocontigs", assembly="nanopore.canu", kpneumo.refs,
                     parsnp.label="nano_parsnp", full.tree=F) {

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
        file.copy(kpneumo.refs$dest.name, file.path(contig.dir, kpneumo.refs$local.name))
        system(paste0("gunzip -f ", file.path(contig.dir, "*.gz")))
        kpneumo.refs$local.name=gsub(".gz", "", kpneumo.refs$local.name)
    } else {
        ##Just HS11286 for now
        file.copy(kpneumo.refs$dest.name[kpneumo.refs$coreref], file.path(contig.dir, kpneumo.refs$local.name[kpneumo.refs$coreref]))
        ##Got to ungzip these for parsnp (obnoxiously)
        system(paste0("gunzip -f ", file.path(contig.dir, "*.gz")))
        kpneumo.refs$local.name=gsub(".gz", "", kpneumo.refs$local.name)
    }
    ##parsnp
    system(paste0("~/Code/Harvest-Linux64-v1.1.2/parsnp -c -r ", file.path(contig.dir, kpneumo.refs$local.name[kpneumo.refs$coreref]),
                  " -d ", contig.dir, " -p 10 -o ./", parsnp.label))
    
    ##Count snps
    vcfloc=file.path(workdir,parsnp.label, "parsnp.vcf")  
    system(paste0("~/Code/Harvest-Linux64-v1.1.2/harvesttools -i ", file.path(workdir, parsnp.label, "parsnp.ggr"), " -V ", vcfloc))
    count.snps(vcfloc)

    ##XMFA
    system(paste0("~/Code/Harvest-Linux64-v1.1.2/harvesttools -i ", file.path(workdir, parsnp.label, "parsnp.ggr"), " -M ",
                  file.path(workdir, parsnp.label, "parsnp.mult.fa")))
    
    ##copy tree file out 
    system(paste0("cp -R ", file.path(workdir, parsnp.label), " ", file.path(outdir)))
}


count.snps = function(vcfloc) {
    ##Count snps in vcf file

    vcfs=read_tsv(vcfloc, comment="##")

    ##ok - for each sample column, need to look at delta for all other sample columns, then sum to get that square on the table

    ##Interesting names are names(vcfs)[10:length(vcfs)]

    ##Sometimes columns are 2s
    comp=as_data_frame(t(combn(names(vcfs[10:length(vcfs)]),2))) %>%
        mutate(val=colSums((vcfs[V1]-vcfs[V2]) != 0 ))

    comp=rbind(comp, data_frame(V1=comp$V2, V2=comp$V1, val=comp$val))

    write_csv(spread(comp, V2, val), paste0(vcfloc, ".csv"))
}

#ref.cols=c("assembly_accession", "bioproject", "biosample", "wgs_master", "refseq_category", "taxid", "species_taxid", "organism_name",
#           "infraspecific_name", "isolate", "version_status", "assembly_level", "release_type", "genome_rep", "seq_rel_date", "asm_name",
#           "submitter", "gbrs_paired_asm", "paired_asm_comp", "ftp_path", "excluded_from_refseq")

#kpneumo.refs=read_tsv(file.path(outdir, "assembly_summary.txt"), comment="#")

##All kpneumo refs
kpneumo.refs=read_csv(file.path("~/Dropbox/Data/Nanopore/170104_reassemble", "genomes_proks.csv")) %>%
    rename(org.name=`#Organism/Name`) %>%
    mutate(rootname=sapply(strsplit(`GenBank FTP`, split="/"), function(x) {x[length(x)]})) %>%
    mutate(source.name=paste0(`GenBank FTP`, "/", rootname, "_genomic.fna.gz")) %>%
    mutate(dest.name=file.path(refdir, paste0(rootname, "_genomic.fna.gz"))) %>%
    filter(Strain!="NA") %>%
    mutate(Strain=gsub(" ", "", Strain)) %>%
    mutate(Strain=gsub("[[:punct:]]", "", Strain)) %>%
    mutate(coreref=(Strain=="HS11286")) %>%
    mutate(local.name=paste0(Strain, ".fasta.gz"))

    


##Download 
if (FALSE) {

    for (i in 1:dim(kpneumo.refs)[1]) {
        ##Get all refs
        download.file(kpneumo.refs$source.name[i], destfile=kpneumo.refs$dest.name[i])
    }
    
}

setwd(workdir)

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1_WT3RQSVGvR97-asIHtIy0WWg49rAFiWWQe9g8BWNiQ/edit?usp=sharing")
dataloc=gs_read(fullsheet, ws="KPneumo0518")


##Raw nanopore assembly tree

#make.tree(workdir, outdir, label="nanocontigs", assembly="nanopore.canu2", kpneumo.refs, parsnp.label="nano_parsnp")
#make.tree(workdir, outdir, label="piloncontigs", assembly="pilon2", kpneumo.refs, parsnp.label="pilon_parsnp")

##Next

##ok - make tables for illumina only, nanopore raw, nanopore pilonx1, nanopore pilon looped
make.tree(dataloc, workdir, outdir, label="spadescontigs", assembly="illumina.spades", kpneumo.refs, parsnp.label="spades_parsnp")
make.tree(dataloc, workdir, outdir, label="nanoporeraw", assembly="nanopore.canu2", kpneumo.refs, parsnp.label="nanopore_raw_parsnp")
make.tree(dataloc, workdir, outdir, label="nanoporeonepilon", assembly="pilon2", kpneumo.refs, parsnp.label="pilon_one_parsnp")
make.tree(dataloc, workdir, outdir, label="pilonloop", assembly="pilon.cor", kpneumo.refs, parsnp.label="pilon_loop_parsnp")



##Generate massive tree with nanopore assemblies and all K. Pnemuo refs
kpneumo.comp=filter(kpneumo.refs, Level=="Complete Genome")

make.tree(dataloc, workdir, outdir, label="pilonloopcomplete", assembly="pilon.cor", kpneumo.comp,
          parsnp.label="pilon_loop_complete_parsnp", full.tree=T)

make.tree(dataloc, workdir, outdir, label="spadescomplete", assembly="illumina.spades", kpneumo.comp,
          parsnp.label="spades_complete_parsnp", full.tree=T)



##Ok - looking just at location of SNPs between strains in pilon loop

pilonloopvcf=read_tsv(file.path(outdir, "pilon_loop_parsnp", "parsnp.vcf"), comment="##") %>%
    rename(chr=`#CHROM`)

##Tree branches are 1,2,6,7 and 9,8,4,12
XDR=paste0("isolate_", c(9,8,4,10,12), ".fasta")
HMV=paste0("isolate_", c(1,2,6,7), ".fasta")


##Do these - make sense given the number of isolates?
HMV.snp=pilonloopvcf %>%
    select(1:9,one_of(HMV)) %>%
    mutate(diffy=rowSums(cbind((.[10]-.[11])!=0, (.[10]-.[12])!=0, (.[10]-.[13])!=0))) %>%
    filter(diffy>0)

XDR.snp=pilonloopvcf %>%
    select(1:9,one_of(XDR)) %>%
    mutate(diffy=rowSums(cbind((.[10]-.[11])!=0, (.[10]-.[12])!=0, (.[10]-.[13])!=0, (.[10]-.[14])!=0))) %>%
    filter(diffy>0)


HMV.gff=HMV.snp %>%
    mutate(source="parsnp", feature="SNP", strand=".", frame=".", attribute=paste0("ALT=", ALT), start=POS, end=POS) %>%
    select(seqname=chr, source=source, feature=feature, start=start, end=end, score=QUAL, strand=strand, frame=frame, attribute=attribute)

XDR.gff=XDR.snp %>%
    mutate(source="parsnp", feature="SNP", strand=".", frame=".", attribute=paste0("ALT=", ALT), start=POS, end=POS) %>%
    select(seqname=chr, source=source, feature=feature, start=start, end=end, score=QUAL, strand=strand, frame=frame, attribute=attribute)


write_tsv(HMV.gff, file.path(outdir, "HMV.gff"), col_names=F)
write_tsv(XDR.gff, file.path(outdir, "XDR_SNP.gff"), col_names=F)


make.tree(dataloc, workdir, outdir, label="spadescontigs", assembly="illumina.spades", kpneumo.refs, parsnp.label="spades_parsnp")

kpneumo.comp=filter(kpneumo.refs, Level=="Complete Genome")

HMV.ref=tibble(dest.name=dataloc$pilon.cor[dataloc$trish.id==4],
               local.name="isolate_1.fasta",
               coreref=T)

XDR.ref=tibble(dest.name=dataloc$pilon.cor[dataloc$trish.id==1],
               local.name="isolate_4.fasta",
               coreref=T)

make.tree(dataloc[dataloc$trish.id %in% c(9,8,4,10,12),], workdir, outdir, label="XDRparsnp", assembly="pilon.cor", XDR.ref,
          parsnp.label="XDRparsnp", full.tree=F)

make.tree(dataloc[dataloc$trish.id %in% c(1,2,6,7),], workdir, outdir, label="HMVparsnp", assembly="pilon.cor", HMV.ref,
          parsnp.label="HMVparsnp", full.tree=F)


