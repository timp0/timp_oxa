library(tidyverse)
library(googlesheets)

workdir="/atium/Data/Nanopore/Analysis/170516_oxapilon/annot/"
outdir="~/Dropbox/Data/Nanopore/170515_oxa/"


##Spreadsheet with sample information, could also be used from a csv
gs_auth(token = "~/googlesheets_token.rds")

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1_WT3RQSVGvR97-asIHtIy0WWg49rAFiWWQe9g8BWNiQ/edit?usp=sharing")
dataloc=gs_read(fullsheet, ws="KPneumo0518") %>%
    filter(Sequenced=="both")

blast.cols=c("query.id", "subject.id", "per.identity", "align.len", "mismatch", "gap.opens", "q.start", "q.end", "s.start", "s.end", "evalue", "bit.score")

card.dir="/atium/Data/Reference/CARD/"

card.tab=read_csv(file.path(card.dir, "aro.csv")) %>%
    mutate(Description=gsub('\n', '', Description))


##Fix card tab - clean to remove \ns
##Fix to make it one hit only for card blast?

if (!dir.exists(workdir)) {
    dir.create(workdir)
}

if (!dir.exists(outdir)) {
    dir.create(outdir)
}

if (!dir.exists(file.path(outdir, "annot"))) {
    dir.create(file.path(outdir, "annot"))
}


setwd(workdir)



##Load CSV with sample info

for (i in 1:dim(dataloc)[1]) {

    
    samp=dataloc$trish.id[i]
    
    ##Run Prokka
    system(paste0("prokka --outdir ", workdir, "/", samp, "_genus --genus Klebsiella --usegenus --prefix ", samp, "_genus ", dataloc$pilon.cor[i]))
    
    file.copy(file.path(workdir, paste0(samp, "_genus"), paste0(samp, "_genus.gff")), file.path(outdir, "annot", paste0(samp, "_genus.gff")))
    
    ##Run RGI (CARD)
    system(paste0("/bin/bash -c ", shQuote(paste0("source activate py27; rgi -t contig -i ", dataloc$pilon.cor[i],
                                                  " -n 8 -o ", samp, ".rgicard"))))

    system(paste0("/bin/bash -c ", shQuote(paste0("source activate py27; rgi_jsontab -i ", samp, ".rgicard.json -o ", samp, ".rgicard"))))
   
    rgi=read_tsv(file.path(workdir, paste0(samp, ".rgicard.txt"))) %>%
        mutate(clabel=gsub('_[0-9]*$', '', ORF_ID)) %>%
        mutate(attrib=paste0("Name=",Best_Hit_ARO, ";gene=", Best_Hit_ARO, ";ARO_category=", Best_Hit_ARO_category, ";cutoff=", CUT_OFF))
    
    rgi.gff=rgi %>%
        mutate(source="rgicard", feature="CDS", frame="0") %>%
        select(seqname=clabel, source=source, feature=feature, start=START,
               end=STOP, score=Best_Identities, strand=ORIENTATION,
               frame=frame, attribute=attrib) 
    
    write_tsv(rgi.gff, file.path(outdir, "annot", paste0(samp, ".rgicard.gff")), col_names=F)
    
    
    ##Run CARD Blast
    system(paste0("blastn -db /atium/Data/Reference/CARD/CARDblast -query ", dataloc$pilon.cor[i],
                  " -evalue 1e-50 -outfmt 6 -best_hit_overhang .1 -best_hit_score_edge 1e-5 >", workdir, "/", samp, ".cardblast.tsv"))
    
    card.res=read_tsv(file.path(workdir, paste0(samp, ".cardblast.tsv")), col_names=blast.cols) %>%
        mutate(aro=sapply(strsplit(subject.id, split="\\|"), function(x) {x[5]})) %>%
        mutate(aro.idx=pmatch(aro, card.tab$Accession, duplicates.ok=TRUE)) %>%
        mutate(aro.name=card.tab$Name[aro.idx]) %>%
        mutate(aro.desc=card.tab$Description[aro.idx])   
    
    card.gff=card.res %>%
        mutate(source="cardblast", feature="CDS", frame="0") %>%
        mutate(strand=ifelse(s.end>s.start, yes="+", no="-")) %>%
        mutate(attrib=paste0("Name=", aro.name, ";gene=", aro.name, ";per.identity=", per.identity, ";bit.score=", bit.score,";evalue=",evalue, ";ARO=",aro, ";Description=",aro.desc)) %>%
        select(seqname=query.id, source=source, feature=feature, start=q.start, end=q.end, score=bit.score, strand=strand, frame=frame, attribute=attrib)

    write_tsv(card.gff, file.path(outdir, "annot", paste0(samp, ".cardblast.gff")), col_names=F)
    
    
    
    ##Look for specific gene mutations
    mutgenedb="/atium/Data/Nanopore/Analysis/161223_kpneumo/genedb/mutdbase"
    system(paste0("blastn -db ", mutgenedb, " -query ", dataloc$pilon.cor[i], " -outfmt 6 >", workdir, samp, ".mutblast.tsv"))
    
    mut.res=read_tsv(file.path(workdir, paste0(samp, ".mutblast.tsv")), col_names=blast.cols)
    
    mut.gff=mut.res %>%
        mutate(source="mutblast", feature="CDS", frame="0") %>%
        mutate(strand=ifelse(s.end>s.start, yes="+", no="-")) %>%
        mutate(attrib=paste0("Name=", subject.id, ";gene=", subject.id, ";per.identity=", per.identity, ";bit.score=", bit.score,";evalue=",evalue)) %>%
        select(seqname=query.id, source=source, feature=feature, start=q.start, end=q.end, score=bit.score, strand=strand, frame=frame, attribute=attrib)
    
    write_tsv(mut.gff, file.path(outdir, "annot", paste0(samp, ".mutblast.gff")), col_names=F)
    
    
    ##look for specific is recog sequences
    isedb="/atium/Data/Nanopore/Analysis/161223_kpneumo/genedb/isecp1blast"
    system(paste0("blastn -task blastn-short -db ", isedb, " -query ", dataloc$pilon.cor[i], " -outfmt 6 >", workdir, samp, ".iseblast.tsv"))
    
    ise.res=read_tsv(file.path(workdir, paste0(samp, ".iseblast.tsv")), col_names=blast.cols)
    
    ise.gff=ise.res %>%
        mutate(source="iseblast", feature="Mobile_element", frame="0") %>%
        mutate(strand=ifelse(s.end>s.start, yes="+", no="-")) %>%
        mutate(attrib=paste0("Name=", subject.id, ";gene=", subject.id, ";per.identity=", per.identity, ";bit.score=", bit.score,";evalue=",evalue)) %>%
        select(seqname=query.id, source=source, feature=feature, start=q.start, end=q.end, score=bit.score, strand=strand, frame=frame, attribute=attrib)
    
    write_tsv(ise.gff, file.path(outdir, "annot", paste0(samp, ".iseblast.gff")), col_names=F)


    ##MLST typing
    system(paste0("mlst ", dataloc$pilon.cor[i], " >", workdir, samp, ".mlst.tsv"))
    
}



