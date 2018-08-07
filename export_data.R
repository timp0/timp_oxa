library(tidyverse)
library(googlesheets)

gs_auth(token = "~/googlesheets_token.rds")

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/10lS7CqnZCWhpsMgaUsZRaWDFgnI_LU8PrnsYRtlgikY/edit?usp=sharing")
dataloc=gs_read(fullsheet, ws="KPneumo180212")

outdir="~/Dropbox/timplab_data/cpowgs/180804_oxa"

gs_auth(token = "~/googlesheets_token.rds")

if (!dir.exists(outdir)) {
    dir.create(outdir)
}

##just 9.4

dataloc=dataloc %>%
    filter(nanopore.run=="R9.4.2D")

file.copy(dataloc$illumina.r1, file.path(outdir, basename(dataloc$illumina.r1)))
file.copy(dataloc$illumina.r2, file.path(outdir, basename(dataloc$illumina.r2)))
file.copy(dataloc$nanopore.fasta, file.path(outdir, basename(dataloc$nanopore.fasta)))
