library(ggplot2)
d = read.table(snakemake@input, header=F)
colnames(d) = c("seqName", "read", "len", "firstBase", "lastBase", "trimmed")
d$readNr = matrix(unlist(strsplit(as.character(d$read), ":")), ncol=4, byrow=T)[,1]

ggplot(d) + geom_histogram(aes(len, fill=readNr), position="dodge", binwidth=5)
ggsave(snakemake@output)
