library(ggplot2)
library(scales)
library(grid)
library(RColorBrewer)
library(reshape2)

args=commandArgs(trailingOnly=TRUE)

x = read.table(args[1], header=TRUE)
if (is.na(args[2])) quit()

txtFontSize=16
axisFontSize=16
axisTtlFontSize=22
lgdTtlFontSize=22
lgdFontSize=16
scienceTheme=theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_blank(), legend.background=element_blank(), panel.background = element_blank(), panel.border=element_blank(), strip.background = element_blank(), axis.line=element_line(size=0.7, color="black"), axis.text.x=element_text(size=axisFontSize), axis.text.y=element_text(size=axisFontSize), axis.title.x=element_text(size=axisTtlFontSize), axis.title.y=element_text(size=axisTtlFontSize), legend.title=element_text(size=lgdTtlFontSize, face="bold"), legend.text=element_text(size=lgdFontSize), text=element_text(size=txtFontSize))

png(paste0(args[2], ".png"), width=800, height=400)
subset=x[x$id==args[2],]
chrName = as.character(unique(x$chr)[1])
p1=ggplot(data=subset, aes(x=start, y=wratio)) + geom_line(aes(color=type), size=0.1) + geom_point(aes(color=type, size=support))
p1=p1 + scienceTheme + xlab(chrName) + ylab("Watson Ratio") 
p1=p1 + scale_x_continuous(labels=comma) + ylim(0,1) + theme(axis.text.x=element_text(angle=45, hjust=1))
p1
dev.off()
print(warnings())
