library(tidyverse)
library(googlesheets)

workdir="/atium/Data/Nanopore/Analysis/180114_oxa"
outdir="~/Dropbox/Data/Nanopore/180114_oxa"
refdir=file.path(workdir, "ref")
dir.create(refdir)


gs_auth(token = "~/googlesheets_token.rds")

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1_WT3RQSVGvR97-asIHtIy0WWg49rAFiWWQe9g8BWNiQ/edit?usp=sharing")
dataloc=gs_read(fullsheet, ws="KPneumo0518")

setwd(workdir)

##From big tree code - closest refs are:
dataloc$close.ref=NA
dataloc$close.ref[dataloc$trish.id %in% c(1, 2, 7, 6)]="ED23"
dataloc$close.ref[dataloc$trish.id %in% c(4,8,9,10,12)]="MS6671"

##All kpneumo refs
kpneumo.refs=read_csv(file.path("~/Dropbox/Data/Nanopore/170104_reassemble", "genomes_proks.csv")) %>%
    rename(org.name=`#Organism/Name`) %>%
    mutate(rootname=sapply(strsplit(`GenBank FTP`, split="/"), function(x) {x[length(x)]})) %>%
    mutate(source.name=paste0(`GenBank FTP`, "/", rootname, "_genomic.fna.gz")) %>%
    filter(Strain!="NA") %>%
    mutate(Strain=gsub(" ", "", Strain)) %>%
    mutate(Strain=gsub("[[:punct:]]", "", Strain)) %>%
    filter(Strain %in% c("ED23", "HS11286", "MS6671")) %>%
    mutate(local.name=paste0(Strain, ".fasta"))

##Download 
if (FALSE) {
    for (i in 1:dim(kpneumo.refs)[1]) {
        ##Get all refs
        download.file(kpneumo.refs$source.name[i], destfile=file.path(refdir, paste0(kpneumo.refs$local.name[i], ".gz")))
        system(paste0("gunzip ",file.path(refdir, paste0(kpneumo.refs$local.name[i], ".gz"))))
    }    
}

##Index

if (FALSE) {

    for (i in 1:dim(kpneumo.refs)[1]) {
        ##Index with bowtie2
        system(paste0("bowtie2-build ", file.path(refdir, kpneumo.refs$local.name[i]), " ", file.path(refdir, kpneumo.refs$Strain[i])))
        ##Index with bwa
        system(paste0("bwa index ", file.path(refdir, kpneumo.refs$local.name[i])))
    }

}

dataloc=dataloc %>%
    mutate(final.btbam=ifelse(is.na(pilon.cor), "NA", file.path(workdir, paste0(trish.id, ".sorted.bt2.bam")))) %>%
    mutate(final.bwabam=ifelse(is.na(pilon.cor), "NA", file.path(workdir, paste0(trish.id, ".sorted.bwa.bam")))) 


##Align
if (FALSE) {
    for (i in 1:dim(dataloc)[1]) {
        if (!is.na(dataloc$close.ref[i])) {
            ##Align ill to final bowtie2
            this.ref=which(kpneumo.refs$Strain==dataloc$close.ref[i])
            system(paste0("bowtie2 -p 15 -x ", file.path(refdir, kpneumo.refs$Strain[this.ref]), " -1 ", dataloc$illumina.r1[i], " -2 ",
                          dataloc$illumina.r2[i], " | samtools view -bS - | samtools sort - -o ", dataloc$final.btbam[i]))
            system(paste0("samtools index ", dataloc$final.btbam[i]))
            
            ##Align with bwa nanopore to final
            if (dataloc$Sequenced[i] == "both") { 
                system(paste0("~/Code/bwa/bwa mem -t 15 -x ont2d ", file.path(refdir, kpneumo.refs$local.name[this.ref]), " ", dataloc$nanopore.fasta[i],
                              " | samtools view -S -b - | samtools sort - -o ", dataloc$final.bwabam[i]))
                system(paste0("samtools index ", dataloc$final.bwabam[i]))
            }
        }
    }
}


##Ok - get consensus from bam(s)


dataloc=dataloc %>%
    mutate(final.ill.cons=ifelse(is.na(pilon.cor), "NA", file.path(workdir, paste0(trish.id, ".ill")))) %>%
    mutate(final.nano.cons=ifelse(is.na(pilon.cor), "NA", file.path(workdir, paste0(trish.id, ".nano")))) 


if (FALSE) {
    for (i in 1:dim(dataloc)[1]) {
        if (!is.na(dataloc$close.ref[i])) {
            this.ref=which(kpneumo.refs$Strain==dataloc$close.ref[i])
            system(paste0("samtools mpileup -uf ", file.path(refdir, kpneumo.refs$local.name[this.ref]),
                          " ", dataloc$final.btbam[i], " | bcftools call -c - | vcfutils.pl vcf2fq > ",
                          dataloc$final.ill.cons[i], ".fq"), wait=F)


            system(paste0("samtools mpileup -uf ", file.path(refdir, kpneumo.refs$local.name[this.ref]),
                          " ", dataloc$final.bwabam[i], " | bcftools call -c - | vcfutils.pl vcf2fq > ",
                          dataloc$final.nano.cons[i], ".fq"), wait=F)

        }
    }
}


fullsheet=gs_url("https://docs.google.com/spreadsheets/d/10lS7CqnZCWhpsMgaUsZRaWDFgnI_LU8PrnsYRtlgikY/edit?usp=sharing")
gs_ws_new(fullsheet, ws_title="KPneumo180212", input=dataloc)

for (i in 1:dim(dataloc)[1]) {
    if (!is.na(dataloc$close.ref[i])) {
        
        system(paste0("seqtk seq -a ", dataloc$final.ill.cons[i], ".fq >", dataloc$final.ill.cons[i], ".fa"))
        system(paste0("seqtk seq -a ", dataloc$final.nano.cons[i], ".fq >", dataloc$final.nano.cons[i], ".fa"))

    }
}
