library(tidyverse)
library(googlesheets)

workdir="/atium/Data/Nanopore/Analysis/170330_oxa"
outdir="~/Dropbox/Data/Nanopore/170515_oxa"

gs_auth(token = "~/googlesheets_token.rds")

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1_WT3RQSVGvR97-asIHtIy0WWg49rAFiWWQe9g8BWNiQ/edit?usp=sharing")
dataloc=gs_read(fullsheet, ws="KPneumo0402")

setwd(workdir)

dir.create(file.path(workdir, "pilonidx"))

dataloc=dataloc %>%
    mutate(pilon.idxfa=ifelse(is.na(pilon2), "NA", file.path(workdir, "pilonidx", paste0(trish.id, ".pilon.fa")))) %>%
    mutate(pilon.btidx=ifelse(is.na(pilon2), "NA", file.path(workdir, "pilonidx", paste0(trish.id, ".pilon")))) %>%
    mutate(pilon.bwaidx=ifelse(is.na(pilon2), "NA", file.path(workdir, "pilonidx", paste0(trish.id, ".pilon"))))



if (TRUE) {
    for (i in 1:dim(dataloc)[1]) {

        if (!is.na(dataloc$pilon.idxfa[i])) {
            ##Index canu with bowtie2

            file.copy(dataloc$pilon2[i], dataloc$pilon.idxfa[i])
            
            system(paste0("bowtie2-build ", dataloc$pilon.idxfa[i], " ", dataloc$pilon.btidx[i]))
            system(paste0("~/Code/bwa/bwa index ", dataloc$pilon.idxfa[i]))
        }
    }   
}



##align illumina to nanopore contigs
dir.create(file.path(workdir, "pilonbam"))

dataloc=dataloc %>%
    mutate(pilon.btbam=ifelse(is.na(pilon2), "NA", file.path(workdir, "pilonbam", paste0(trish.id, ".sorted.bt2.bam")))) %>%
    mutate(pilon.bwabam=ifelse(is.na(pilon2), "NA", file.path(workdir, "pilonbam", paste0(trish.id, ".sorted.bwa.bam")))) 

        

if (TRUE) {

    for (i in 1:dim(dataloc)[1]) {
        if (!is.na(dataloc$pilon.idxfa[i])) {
            ##Align ill to pilon bowtie2
            system(paste0("bowtie2 -p 15 -x ", dataloc$pilon.btidx[i], " -1 ", dataloc$illumina.r1[i], " -2 ",
                          dataloc$illumina.r2[i], " | samtools view -bS - | samtools sort - -o ", dataloc$pilon.btbam[i]))
            system(paste0("samtools index ", dataloc$pilon.btbam[i]))

            ##Align with bwa nanopore to pilon
            system(paste0("~/Code/bwa/bwa mem -t 15 -x ont2d ", dataloc$pilon.idxfa[i], " ", dataloc$nanopore.fasta[i], " | samtools view -S -b - | samtools sort - -o ",
                          dataloc$pilon.bwabam[i]))
            system(paste0("samtools index ", dataloc$pilon.bwabam[i]))
            

        }
    }
}


if (TRUE) {
    vcfloc="/home/timp/Dropbox/Data/Nanopore/170402_reassemble/pilon_parsnp/parsnp.vcf"

    vcfs=read_tsv(vcfloc, comment="##")




#fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1_WT3RQSVGvR97-asIHtIy0WWg49rAFiWWQe9g8BWNiQ/edit?usp=sharing")
#gs_ws_new(fullsheet, ws_title="KPneumo0402", input=dataloc)
