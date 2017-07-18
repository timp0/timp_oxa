library(tidyverse)
library(googlesheets)

workdir="/atium/Data/Nanopore/Analysis/170516_oxapilon"
outdir="~/Dropbox/Data/Nanopore/170515_oxa"

gs_auth(token = "~/googlesheets_token.rds")

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1_WT3RQSVGvR97-asIHtIy0WWg49rAFiWWQe9g8BWNiQ/edit?usp=sharing")
dataloc=gs_read(fullsheet, ws="KPneumo0518")

setwd(workdir)

dir.create(file.path(workdir, "assembly.coverage"))
dir.create(file.path(workdir, "reference.coverage"))
dataloc=dataloc %>%
    mutate(final.idxfa=ifelse(is.na(pilon.cor), "NA", file.path(workdir, "assembly.coverage", paste0(trish.id, ".loop.fa")))) %>%
    mutate(final.btidx=ifelse(is.na(pilon.cor), "NA", file.path(workdir, "assembly.coverage", paste0(trish.id, ".loop")))) %>%
    mutate(final.bwaidx=ifelse(is.na(pilon.cor), "NA", file.path(workdir, "asssembly.coverage", paste0(trish.id, ".loop"))))
    
ref.bwaidx="/mithril/Data/NGS/Reference/kpneumo/NC_106845.fasta"
ref.btidx="/mithril/Data/NGS/Reference/kpneumo/kpneumo"


if (TRUE) {
    for (i in 1:dim(dataloc)[1]) {

        if (!is.na(dataloc$final.idxfa[i])) {
            ##Index canu with bowtie2

            file.copy(dataloc$pilon.cor[i], dataloc$final.idxfa[i])
            
            system(paste0("bowtie2-build ", dataloc$final.idxfa[i], " ", dataloc$final.btidx[i]))
            system(paste0("~/Code/bwa/bwa index ", dataloc$final.idxfa[i]))
        }
    }   
}



##align illumina to nanopore contigs
dataloc=dataloc %>%
    mutate(final.btbam=ifelse(is.na(pilon.cor), "NA", file.path(workdir, "assembly.coverage", paste0(trish.id, ".sorted.bt2.bam")))) %>%
    mutate(final.bwabam=ifelse(is.na(pilon.cor), "NA", file.path(workdir, "assembly.coverage", paste0(trish.id, ".sorted.bwa.bam")))) 

        

if (TRUE) {

    for (i in 1:dim(dataloc)[1]) {
        if (!is.na(dataloc$final.idxfa[i])) {
            ##Align ill to final bowtie2
            system(paste0("bowtie2 -p 15 -x ", dataloc$final.btidx[i], " -1 ", dataloc$illumina.r1[i], " -2 ",
                          dataloc$illumina.r2[i], " | samtools view -bS - | samtools sort - -o ", dataloc$final.btbam[i]))
            system(paste0("samtools index ", dataloc$final.btbam[i]))

            ##Align with bwa nanopore to final
            system(paste0("~/Code/bwa/bwa mem -t 15 -x ont2d ", dataloc$final.idxfa[i], " ", dataloc$nanopore.fasta[i], " | samtools view -S -b - | samtools sort - -o ",
                          dataloc$final.bwabam[i]))
            system(paste0("samtools index ", dataloc$final.bwabam[i]))
            

        }
    }
}

