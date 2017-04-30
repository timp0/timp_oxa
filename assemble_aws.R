library(tidyverse)

##Run canu
if (TRUE) {

    dataloc=tibble(trish.id=c(1,2,4,6,7,8,9,10,12), nanopore.fasta=file.path("nanopore.raw", paste0(c(1,2,4,6,7,8,9,10,12), ".nanopore.2D.fasta")))
    
    dataloc=dataloc %>%
        mutate(nanopore.canu=NA)
    for (i in 1:dim(dataloc)[1]) {

        samp=dataloc$trish.id[i]

        ##Original assembly
        ##system(paste0("~/canu/Linux-amd64/bin/canu -p ", samp, ".nanopore -d ", paste0("canu/", samp), ".nanopore -genomeSize=5m -nanopore-raw ", dataloc$nanopore.fasta[i]))
        ##plasmid assembly
        system(paste0("~/canu/Linux-amd64/bin/canu -p ", samp, ".nanopore -d ", paste0("canu/", samp), ".nanopore -genomeSize=5m -contigFilter=\"2 1000 1.0 1.0 2\" -corOutCoverage=1000 -nanopore-raw ", dataloc$nanopore.fasta[i]))
              
        #dataloc$nanopore.canu[i]=file.path(workdir, paste0(samp, ".nanopore"), paste0(samp, ".nanopore.contigs.fasta"))
        
    }

}

if (FALSE) {

    dataloc=tibble(trish.id=1:12) %>%
        mutate(illumina.r1=file.path("illumina.raw", paste0(trish.id, "_R1.fastq.gz"))) %>%
        mutate(illumina.r2=file.path("illumina.raw", paste0(trish.id, "_R2.fastq.gz")))

    for (i in 1:dim(dataloc)[1]) {

        samp=dataloc$trish.id[i]

        system(paste0("~/SPAdes-3.10.0/bin/spades.py -t 30 -1 ", dataloc$illumina.r1[i],
                      " -2 ", dataloc$illumina.r2[i], " -o ", "spades/", samp))

        #dataloc$illumina.spades[i]=file.path(workdir, paste0(samp, ".spades"))
    }        

    
}
