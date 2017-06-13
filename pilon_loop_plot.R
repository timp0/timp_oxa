library(tidyverse)

##wc -l *changes >pilon.changes.txt

pilondir="/atium/Data/Nanopore/Analysis/170516_oxapilon/pilon"
plotdir="~/Dropbox/Data/Nanopore/170515_oxa"


meta=read_table(file.path(pilondir, "pilon.changes.txt"), col_names=F) %>%
    filter(X2 != "total")

colnames(meta)=c("num", "fname")

meta=meta %>%
    separate(fname, c("samp", "round"), remove=F, extra="drop", convert=T)

pdf(file.path(plotdir, "pilons.pdf"))

ggplot(meta, aes(x=round, y=num, group=factor(samp), color=factor(samp)))+theme_classic()+geom_line()+scale_y_log10()

dev.off()
