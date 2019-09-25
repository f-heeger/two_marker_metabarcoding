library(ggplot2)
        
raw=read.table(snakemake@input[[1]])
colnames(raw) = c("sample", "x", "y")

allPlotData = data.frame()
anno=data.frame()

for (tSample in unique(raw$sample)){

    print(tSample)
    d = subset(raw, sample==tSample)
    plotData = data.frame(unique(d$x))
    colnames(plotData) = c("x")
    N=length(plotData$x)

    plotData$mean=NA
    plotData$cmin=NA
    plotData$cmax=NA
    plotData$sample=tSample

    conf = 0.95

    for (i in 1:length(plotData$x)) {
        x_i = plotData$x[i]
        s = subset(d, x==x_i)
        plotData$mean[i] = mean(s$y)
        plotData$cmin[i] = sort(s$y)[floor(length(s$y)*(1-conf)/2)]
        plotData$cmax[i] = sort(s$y)[ceiling(length(s$y)*(1-(1-conf)/2))]
    }
    anno=rbind(anno, plotData[plotData$x==max(plotData$x),])

allPlotData = rbind(allPlotData, plotData)

}

ggplot(allPlotData, aes(x, mean)) + geom_point(aes(color=sample), size=0.1) + geom_ribbon(aes(ymin=cmin, ymax=cmax, fill=sample), alpha=0.2) + annotate("text", x=anno$x, y=anno$mean, label=anno$sample)
ggsave(snakemake@output[[1]], width=16, height=9)
