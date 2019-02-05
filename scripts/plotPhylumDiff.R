library(reshape2)
library(ggplot2)

d = read.table(snakemake@input, header=T)

newL = unique(c(levels(d$phylum58s), levels(d$phylumIts2)))
d$phylum58s = factor(d$phylum58s, levels=newL)
d$phylumIts2 = factor(d$phylumIts2, levels=newL)

m = melt(d[d$size>1 & d$phylumIts2 == "None",], id.vars=c("otu", "size"))

cPalette = c("grey", rev(c("#cab2d6", "#fdbf6f", "#fb9a99", "#b2df8a", "#a6cee3", "#6a3d9a", "#ff7f00", "#e31a1c", "#33a02c", "#1f78b4", "#b15928")))

ggplot(m[m$variable=="phylum58s",]) + geom_bar(aes(1,fill=value)) + scale_fill_manual(values=cPalette) + theme_bw()
ggsave(snakemake@output, width=5, height=8)
