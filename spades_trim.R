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

if (!dir.exists(outdir)) {
    dir.create(outdir)
}


setwd(workdir)

for (i in 1:1) {
##Fix here
    system(paste0("atropos --aligner insert -T 6 -a CTGTCTCTTATACACATCT ",
                  "-A CTGTCTCTTATACACATCT ",
                  "-q 20 -o /tmp/read1.trim.fq.gz -p /tmp/read2.trim.fq.gz ",
                  "-pe1 ",dataloc$illumina.r1[i]," -pe2 ", dataloc$illumina.r2[i]))
    
    system(paste0("spades.py -1 /tmp/read1.trim.fq.gz -2 /tmp/read2.trim.fq.gz ",
                  " -o ", dataloc$trish.id[i], "/spades"))
}
