library(tidyverse)
library(googlesheets)

gs_auth(token = "~/googlesheets_token.rds")

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/10lS7CqnZCWhpsMgaUsZRaWDFgnI_LU8PrnsYRtlgikY/edit?usp=sharing")
dataloc=gs_read(fullsheet, ws="KPneumo180212")

workdir="/atium/Data/Nanopore/Analysis/180703_oxa"
outdir="~/Dropbox/timplab_data/cpowgs/180703_oxa"

gs_auth(token = "~/googlesheets_token.rds")

if (!dir.exists(workdir)) {
    dir.create(workdir)
}


setwd(workdir)

for (i in c(3, 5, 11)) {
    ##Fix here
    ##tried TCGTCGGCAGCGTC but still contaminated.  Trying this other one now
    file.remove("/tmp/read1.trim.fq.gz")
    file.remove("/tmp/read2.trim.fq.gz")
    #system(paste0("atropos -T 6 -a GTCTCGTGGGCTCGG ",
    #              "-A GTCTCGTGGGCTCGG ",
    #              "-q 20 -o /tmp/read1.trim.fq.gz -p /tmp/read2.trim.fq.gz ",
    #              "-pe1 ",dataloc$illumina.r1[i]," -pe2 ", dataloc$illumina.r2[i]))

    system(paste0("java -jar ~/Code/Trimmomatic-0.38/trimmomatic-0.38.jar PE -phred33 -threads 6 ",
                  dataloc$illumina.r1[i], " ", dataloc$illumina.r2[i],
                  " /tmp/read1.trim.fq.gz /tmp/read1.unpaired.fq.gz",
                  " /tmp/read2.trim.fq.gz /tmp/read2.unpaired.fq.gz",
                  " ILLUMINACLIP:/home/timp/Code/Trimmomatic-0.38/adapters/NexteraPE-PE.fa:2:30:10",
                  " LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40"))
                  

    
    system(paste0("spades.py -1 /tmp/read1.trim.fq.gz -2 /tmp/read2.trim.fq.gz ",
                  " -o ", dataloc$trish.id[i], "/spades"))
}
 
