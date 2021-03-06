library(tidyverse)
library(googlesheets)

workdir="/atium/Data/Nanopore/Analysis/170104_reassemble"
workdir2="/atium/Data/Nanopore/Analysis/170330_oxa"
outdir="~/Dropbox/Data/Nanopore/170104_reassemble"

gs_auth(token = "~/googlesheets_token.rds")

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1_WT3RQSVGvR97-asIHtIy0WWg49rAFiWWQe9g8BWNiQ/edit?usp=sharing")
dataloc=gs_read(fullsheet, ws=1)

setwd(workdir)

##Load CSV with sample info

#dataloc=read_csv(file.path(outdir, "dataloc.csv")) %>%
#    filter(Sequenced=="both")

dataloc=dataloc %>%
    mutate(nanopore.fasta=ifelse(is.na(nanopore.raw), "NA", file.path(workdir, "nanopore.raw", paste0(trish.id, ".nanopore.2D.fasta"))))


if (FALSE) {   
    for (i in 1:dim(dataloc)[1]) {
    
        if (!is.na(dataloc$nanopore.raw[i])) {
            samp=dataloc$trish.id[i]
            
            ##Run Poretools            
            ##Are there two flowcells?
            two.flow=(!is.na(dataloc$nanopore.2nd[i]))
            
            system(paste0("poretools fasta --type 2D ", dataloc$nanopore.raw[i], " >", dataloc$nanopore.fasta[i]))
            
            if (two.flow) {
                system(paste0("poretools fasta --type 2D ", dataloc$nanopore.2nd[i], " >>", dataloc$nanopore.fasta[i]))
            }
        }
    }
}

##Plot/aquire stats on nanopore fasta/q
if (FALSE) {

    library(Biostrings)

    dataloc=dataloc %>%
        mutate(nanopore.fasta=ifelse(is.na(nanopore.raw), "NA", file.path(workdir, "nanopore.raw", paste0(trish.id, ".nanopore.2D.fasta"))))

    dataloc=dataloc %>%
        mutate(nanopore.medlength=0) %>%
        mutate(nanopore.yield=0)
    
    
    for (i in 1:dim(dataloc)[1]) {
        if (nanopore.fasta != "NA") {
            fa=readDNAStringSet(dataloc$nanopore.fasta[i])
            dataloc$nanopore.medlength=median(width(fa))
            dataloc$nanopore.yield=sum(as.numeric(width(fa)))/1e6
        }     
    }
}


##Run SPAdes
if (FALSE) {
    ##copy illumina data to temp dir for transfer up to AWS
    
    for (i in 1:dim(dataloc)[1]) {
        
        samp=dataloc$trish.id[i]
        
        file.copy(dataloc$illumina.r1[i], file.path(workdir, "illumina.raw", paste0(samp, "_R1.fastq.gz")))
        file.copy(dataloc$illumina.r2[i], file.path(workdir, "illumina.raw", paste0(samp, "_R2.fastq.gz")))
    
    }

}

##Set spades location
dataloc=dataloc %>%
    mutate(illumina.spades=file.path(workdir, "spades", trish.id, "contigs.fasta"))

##Ran on AWS, but to run locally:
if (FALSE) {
    
    
    for (i in 1:dim(dataloc)[1]) {
        
        if (!is.na(dataloc$illumina.r1[i])) {
            
            samp=dataloc$trish.id[i]
            
            system(paste0("~/Code/SPAdes-3.9.1/bin/spades.py -1 ", dataloc$illumina.r1[i], " -2 ", dataloc$illumina.r2[i],
                          " -o ", samp, "/spades"))
            
            
            dataloc$illumina.spades[i]=file.path(workdir, paste0(samp, ".spades"))
            
        }
    }
}



##canu
##set canu directories
##2nd canu run is 

    
dataloc=dataloc %>%
    mutate(nanopore.canu=ifelse(is.na(nanopore.raw), "NA", file.path(workdir, "canu", paste0(trish.id, ".nanopore"),
                                                                     paste0(trish.id, ".nanopore.contigs.fasta")))) %>%
    mutate(nanopore.canu2=ifelse(is.na(nanopore.raw), "NA", file.path(workdir2, paste0(trish.id, ".nanopore.contigs.fasta"))))


    
##Run on aws instead
if (FALSE) {

    for (i in 1:dim(dataloc)[1]) {

        if (!is.na(dataloc$nanopore.raw[i])) {
            
            samp=dataloc$trish.id[i]
            
            system(paste0("~/Code/canu/Linux-amd64/bin/canu -p ", samp, ".nanopore -d ", samp, ".nanopore -genomeSize=5m -nanopore-raw ", dataloc$nanopore.fasta[i]))
        }
    }
}


##make index for alignment to canu contigs
dir.create(file.path(workdir, "btidx"))
dir.create(file.path(workdir2, "btidx"))

dataloc=dataloc %>%
    mutate(canu.btidx=ifelse(is.na(nanopore.raw), "NA", file.path(workdir, "btidx", paste0(trish.id, ".nanopore")))) %>%
    mutate(canu2.btidx=ifelse(is.na(nanopore.raw), "NA", file.path(workdir2, "btidx", paste0(trish.id, ".nanopore"))))

if (FALSE) {
    for (i in 1:dim(dataloc)[1]) {

        if (!is.na(dataloc$nanopore.raw[i])) {
            ##Index canu with bowtie2
            
            system(paste0("bowtie2-build ", dataloc$nanopore.canu[i], " ", dataloc$canu.btidx[i]))

        }
    }   
}

if (TRUE) {
    for (i in 1:dim(dataloc)[1]) {

        if (!is.na(dataloc$nanopore.raw[i])) {
            ##Index canu with bowtie2
            
            system(paste0("bowtie2-build ", dataloc$nanopore.canu2[i], " ", dataloc$canu2.btidx[i]))

        }
    }   
}



##align illumina to nanopore contigs
dir.create(file.path(workdir, "btbam"))
dir.create(file.path(workdir2, "btbam"))

dataloc=dataloc %>%
    mutate(canu.ill.align=ifelse(is.na(nanopore.raw), "NA", file.path(workdir, "btbam", paste0(trish.id, ".sorted.bam")))) %>%
    mutate(canu2.ill.align=ifelse(is.na(nanopore.raw), "NA", file.path(workdir2, "btbam", paste0(trish.id, ".sorted.bam")))) 


if (TRUE) {

    for (i in 1:dim(dataloc)[1]) {
        if (!is.na(dataloc$nanopore.raw[i])) {
            ##Align ill to canu bowtie2
            system(paste0("bowtie2 -p 15 -x ", dataloc$canu.btidx[i], " -1 ", dataloc$illumina.r1[i], " -2 ",
                          dataloc$illumina.r2[i], " | samtools view -bS - | samtools sort - -o ", dataloc$canu.ill.align[i]))
            system(paste0("samtools index ", dataloc$canu.ill.align[i]))

        }
    }
}

if (TRUE) {

    for (i in 1:dim(dataloc)[1]) {
        if (!is.na(dataloc$nanopore.raw[i])) {
            ##Align ill to canu bowtie2
            system(paste0("bowtie2 -p 15 -x ", dataloc$canu2.btidx[i], " -1 ", dataloc$illumina.r1[i], " -2 ",
                          dataloc$illumina.r2[i], " | samtools view -bS - | samtools sort - -o ", dataloc$canu2.ill.align[i]))
            system(paste0("samtools index ", dataloc$canu2.ill.align[i]))

        }
    }
}



##Use pilon to fix nanopore contigs with alginments
dir.create(file.path(workdir, "pilon"))
dir.create(file.path(workdir2, "pilon"))

dataloc=dataloc %>%
    mutate(pilon=ifelse(is.na(nanopore.raw), "NA", file.path(workdir, "pilon", paste0(trish.id, ".pilon.fasta")))) %>%
    mutate(pilon2=ifelse(is.na(nanopore.raw), "NA", file.path(workdir2, "pilon", paste0(trish.id, ".pilon.fasta"))))

##run pilon

setwd(workdir)
if (TRUE) {
    for (i in 1:dim(dataloc)[1]) {
        if (!is.na(dataloc$nanopore.raw[i])) {
            system(paste0("java -Xmx30G -jar ~/Code/pilon/pilon-1.21.jar --genome ", dataloc$nanopore.canu[i],
                          " --frags ", dataloc$canu.ill.align[i], " --output ", dataloc$trish.id[i], ".pilon --outdir pilon"))

        }
    }
}

setwd(workdir2)
if (TRUE) {
    for (i in 1:dim(dataloc)[1]) {
        if (!is.na(dataloc$nanopore.raw[i])) {
            system(paste0("java -Xmx30G -jar ~/Code/pilon/pilon-1.22.jar --genome ", dataloc$nanopore.canu2[i],
                          " --frags ", dataloc$canu2.ill.align[i], " --output ", dataloc$trish.id[i], ".pilon --outdir pilon"))

        }
    }
}



fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1_WT3RQSVGvR97-asIHtIy0WWg49rAFiWWQe9g8BWNiQ/edit?usp=sharing")
gs_ws_new(fullsheet, ws_title="KPneumo0402", input=dataloc)
