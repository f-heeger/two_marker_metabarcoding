library(ggplot2)
raw = read.table(snakemake@input[["raw"]], header=F)
raw = raw[order(raw$V1),]
raw$stage = "raw"
d = raw
primer = read.table(snakemake@input[["primer"]], header=F)
primer = primer[order(primer$V1),]
primer$stage = "primerFound"
d = rbind(d, primer)
trimmed = read.table(snakemake@input[["trimmed"]], header=F)
trimmed = trimmed[order(trimmed$V1),]
trimmed$stage = "trimmed"
d=rbind(d, trimmed)
merged = read.table(snakemake@input[["merged"]], header=F)
merged = merged[order(merged$V1),]
merged$stage = "merged"
d=rbind(d,merged)
colnames(d) = c("sample", "readNum", "stage")
d$stage = factor(d$stage, levels=c("merged", "trimmed", "primerFound", "raw"))
ggplot(d) + geom_bar(aes(sample, readNum, fill=stage), stat="identity", position="dodge") + coord_flip() + geom_hline(yintercept = 10000, linetype="dashed")
ggsave(snakemake@output[[1]], width=7, height=28, units="in")
