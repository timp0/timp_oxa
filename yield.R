library(tidyverse)
library(googlesheets)

gs_auth(token = "~/googlesheets_token.rds")

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/10lS7CqnZCWhpsMgaUsZRaWDFgnI_LU8PrnsYRtlgikY/edit?usp=sharing")
dataloc=gs_read(fullsheet, ws="KPneumo180212")

dataloc$ill.r1.yield=0
dataloc$ill.r2.yield=0
dataloc$nano.yield=0


for (i in 1:dim(dataloc)[1]) {
          
    dataloc$ill.r1.yield[i]=as.numeric(system(paste0("gunzip -c ", dataloc$illumina.r1[i],
                                       " | awk 'BEGIN { t=0.0; } ; NR%4==2 {t+=length($0);}END{printf(\"%d\", t);}'"),
                                intern=T))

    dataloc$ill.r2.yield[i]=as.numeric(system(paste0("gunzip -c ", dataloc$illumina.r2[i],
                                       " | awk 'BEGIN { t=0.0; } ; NR%4==2 {t+=length($0);}END{printf(\"%d\", t);}'"),
                                intern=T))

    if (!is.na(dataloc$nanopore.fasta[i])) {
        dataloc$nano.yield[i]=as.numeric(system(paste0("cat ", dataloc$nanopore.fasta[i],
                                            " | awk 'BEGIN { t=0.0; } ; NR%2==0 {t+=length($0);}END{printf(\"%d\", t);}'"),
                                     intern=T))
    }
    
}

gsize=5.33e6

ill.cov=(as.numeric(dataloc$ill.r1.yield)+as.numeric(dataloc$ill.r2.yield))/gsize

nano.cov=as.numeric(dataloc$nano.yield)/gsize


