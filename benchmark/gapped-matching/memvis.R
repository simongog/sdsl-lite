library(ggplot2)
library(reshape2)
library(sitools)
library(scales)

args <- commandArgs(trailingOnly = TRUE)
data <- read.csv(file=args[1],sep=";",header=TRUE);
tdata <- melt(data,id=c("time_ms","event"))
tdata <- subset(tdata,variable!="pid")
dup <-tdata[!duplicated(tdata$event),]

plot <- ggplot(tdata,aes(time_ms,value,color=variable))
plot <- plot + geom_line(size=1)
#plot <- plot + scale_y_log10(expand=c(0,0),labels=f2si,breaks = trans_breaks("log10", function(x) 10^x),"Memory Usage [Byte]")
plot <- plot + scale_y_continuous(expand=c(0,0),labels=f2si,"Memory Usage [Byte]")
#plot <- plot + annotation_logticks(sides = "lr")
plot <- plot + scale_x_continuous(expand=c(0,0),labels=f2si,"Time [Millseconds]")
plot <- plot + scale_linetype(name="Events") + scale_color_discrete(name="Metric")
plot <- plot + geom_vline(data=dup,aes(xintercept = time_ms,linetype=event),show_guide=TRUE)
plot <- plot + theme_grey() + ggtitle("Memory Usage over Time")
ggsave(plot,file=paste0(args[1],".png"))
