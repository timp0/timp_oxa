library(tidyverse)
library(googlesheets)

workdir="/atium/Data/Nanopore/Analysis/180118_oxa"
outdir="~/Dropbox/Data/Nanopore/180114_oxa"
refdir=file.path(workdir, "ref")
dir.create(refdir)


gs_auth(token = "~/googlesheets_token.rds")

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1_WT3RQSVGvR97-asIHtIy0WWg49rAFiWWQe9g8BWNiQ/edit?usp=sharing")
dataloc=gs_read(fullsheet, ws="KPneumo0518")

setwd(workdir)

