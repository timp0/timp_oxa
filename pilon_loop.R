library(tidyverse)
library(googlesheets)

workdir="/atium/Data/Nanopore/Analysis/170516_oxapilon"

outdir="~/Dropbox/Data/Nanopore/170515_oxa"

gs_auth(token = "~/googlesheets_token.rds")

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1_WT3RQSVGvR97-asIHtIy0WWg49rAFiWWQe9g8BWNiQ/edit?usp=sharing")
dataloc=gs_read(fullsheet, ws="KPneumo0402")

setwd(workdir)

##make index for alignment to canu contigs
dir.create(file.path(workdir, "btidx"))

##intially corrected is equal to canu assembly
dataloc=dataloc %>%
    mutate(pilon.cor=ifelse(is.na(nanopore.raw), "NA", nanopore.canu2)) %>%
    mutate(pilon.round=0)


##Just run one do this:
##dataloc=filter(dataloc, trish.id==12)


sttime=proc.time()

for (round in 1:10){

    dataloc=dataloc %>%
        mutate(pilon.btidx=ifelse(is.na(nanopore.raw), "NA", file.path(workdir, "btidx", paste0(trish.id, ".", round, ".nanopore"))))

    for (i in 1:dim(dataloc)[1]) {
        ##for (i in 1:1) {
        if (!is.na(dataloc$nanopore.raw[i])) {
            ##Index canu with bowtie2            
            system(paste0("bowtie2-build ", dataloc$pilon.cor[i], " ", dataloc$pilon.btidx[i]))
            
        }
    }   
    
    
    
    ##align illumina to nanopore contigs
    dir.create(file.path(workdir, "btbam"))
    
    dataloc=dataloc %>%
        mutate(pilon.ill.align=ifelse(is.na(nanopore.raw), "NA", file.path(workdir, "btbam", paste0(trish.id, ".", round, ".sorted.bam"))))
    
    for (i in 1:dim(dataloc)[1]) {
        ##for (i in 1:1) {
        if (!is.na(dataloc$nanopore.raw[i])) {
            ##Align ill to canu bowtie2
            system(paste0("bowtie2 -p 15 -x ", dataloc$pilon.btidx[i], " -1 ", dataloc$illumina.r1[i], " -2 ",
                          dataloc$illumina.r2[i], " | samtools view -bS - | samtools sort -o ", dataloc$pilon.ill.align[i]))
            system(paste0("samtools index ", dataloc$pilon.ill.align[i]))
            
        }
    }

    
    ##Use pilon to fix nanopore contigs with alginments
    dir.create(file.path(workdir, "pilon"))

    oldtigs=dataloc$pilon.cor

    
    dataloc=dataloc %>%
        mutate(pilon.round=round) %>%
        mutate(pilon.cor=ifelse(is.na(nanopore.raw), "NA", file.path(workdir, "pilon", paste0(trish.id,".", pilon.round, ".pilon.fasta"))))
    
    pilon.prefix=paste0(dataloc$trish.id, ".", dataloc$pilon.round, ".pilon")
    
    
    
    ##run pilon
    
    setwd(workdir)
    if (TRUE) {
        for (i in 1:dim(dataloc)[1]) {
        ##for (i in 1:1) {
            if (!is.na(dataloc$nanopore.raw[i])) {
                system(paste0("java -Xmx30G -jar ~/Code/pilon/pilon-1.22.jar --changes --tracks --genome ", oldtigs[i],
                          " --frags ", dataloc$pilon.ill.align[i], " --output ", pilon.prefix[i], " --outdir pilon"))
                
            }
        }
    }
}


for (i in 1:dim(dataloc)[1]) {

    if (!is.na(dataloc$pilon.cor[i])) {
        system(paste0("sed -i -e 's/_pilon//g' ", dataloc$pilon.cor[i]))
    }
}

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1_WT3RQSVGvR97-asIHtIy0WWg49rAFiWWQe9g8BWNiQ/edit?usp=sharing")
gs_ws_new(fullsheet, ws_title="KPneumo0518", input=dataloc)
    
