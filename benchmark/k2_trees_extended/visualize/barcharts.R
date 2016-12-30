require(tikzDevice)
library(gdata)
library(xtable)
library(ggplot2)
library(reshape2)

source("../../basic_functions.R")
scale_data <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/scale_uk2014_morton_num_carray.log")
scale_data$constructs_space <- scale_data$constructs_space/1024/1024/1024
scale_data$constructs_space_vmem <- scale_data$constructs_space_vmem/1024/1024/1024
scale_data$construct_sort_duration <- scale_data$construct_sort_duration/1000
scale_data$construct_duration <- scale_data$construct_duration/1000
scale_data$buildvec_duration <- scale_data$buildvec_duration/1000
scale_data$subtree_constructor_duration <- scale_data$subtree_constructor_duration/1000
scale_data$compression_space <- scale_data$compression_space/1024/1024/1024
scale_data$compression_space_vmem <- scale_data$compression_space_vmem/1024/1024/1024
scale_data$constructor_call_duration <- scale_data$constructor_call_duration/1000
scale_data$comp_word_it <- scale_data$comp_word_it/1000
scale_data$construct_morton_duration <- scale_data$construct_morton_duration/1000
scale_data$frequency_encoding <- scale_data$frequency_encoding/1000
scale_data$dac_compression <- scale_data$dac_compression/1000

scale_2 <- subset(scale_data, select=c('constructor_call_duration', 'threads'))
scale_2 <- scale_2[5:12,]

constructor_call_duration <- melt(scale_2, id.var="threads")
ggplot(constructor_call_duration, aes(x = threads, y = value, fill = variable)) + 
  geom_bar(stat = "identity") + xlab("threads") + ylab("construction time in seconds")+guides(fill=FALSE)

scale_3 <- subset(scale_data, select=c('construct_morton_duration', 'construct_sort_duration', 'construct_duration', 'buildvec_duration', 'threads'))
scale_3 <- scale_3[c(3,4,5,6,7,8,9,11),]
stacked_construction <- melt(scale_3, id.var="threads")
tikz(file = "plot_test.tex", width = 7, height = 4)
plot <- ggplot(stacked_construction, aes(x = threads, y = value, fill = variable, ylim = 10000)) + 
  geom_bar(stat = "identity") + labs(title = "") + xlab("$x$ = threads") + ylab("$y$ = construction time [s]")  +
  scale_x_continuous(breaks = c(1,4,6,8,10,20,30,40))+
  scale_y_continuous(breaks = c(1000,2000,4000,6000,8000,10000,12000,14000))+
  theme(axis.title.y = element_text( vjust=-1.0 ))+
  scale_fill_manual(name = '',
                    labels = c('Morton Number', 
                               'Sort',
                               'Construction',
                               "Vector Merge"),values = c("#E69F00","#D55E00", "#0072B2", "#009E73"))
print(plot)
dev.off()
library(plotrix)
gap.barplot(as.matrix(scale_3),gap=c(7500,12000),xtics=c(0,10,20,30,40))

scale_4 <- subset(scale_data, select=c('comp_word_it', 'frequency_encoding', 'dac_compression', 'threads'))
scale_4 <- scale_4[5:12,]
compression_stacked <- melt(scale_4, id.var="threads")
ggplot(compression_stacked, aes(x = threads, y = value, fill = variable)) + 
  geom_bar(stat = "identity")

scale_5 <- subset(scale_data, select=c('compression_time', 'threads'))
scale_5 <- scale_5[5:12,]
compression <- melt(scale_5, id.var="threads")
ggplot(compression, aes(x = threads, y = value, fill = variable)) + 
  geom_bar(stat = "identity") + xlab("threads") + ylab("compression time in seconds")+guides(fill=FALSE)


scale_6 <- subset(scale_data, select=c('construct_morton_duration', 'construct_sort_duration', 'construct_duration', 'buildvec_duration', 'constructor_call_duration', 'threads'))
scale_6 <- scale_6[3:11,]
one_thread <- subset(scale_6, threads <= 1, )
scale_6$construct_morton_duration <-  one_thread$construct_morton_duration / scale_6$construct_morton_duration
scale_6$construct_sort_duration <-  one_thread$construct_sort_duration / scale_6$construct_sort_duration
scale_6$construct_duration <- one_thread$construct_duration / scale_6$construct_duration
scale_6$buildvec_duration <- one_thread$buildvec_duration / scale_6$buildvec_duration
scale_6$constructor_call_duration <- one_thread$constructor_call_duration / scale_6$constructor_call_duration
scaleab <- melt(scale_6, id.var="threads")

tikz(file = "zord_scale_plot.tex", width = 5, height = 3)
ggplot(data=scaleab, aes(x=threads, y=value, colour=variable)) +
  ylab("$y$ = speedup")+
  xlab("$x$ = threads")+
  geom_line() +
  geom_point()+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_x_continuous(breaks = c(1,10,20,30,40))+
  scale_colour_manual(name = '',
                  labels = c('Morton Number', 
                            'Sort',
                            'Construction',
                              "Vector Merge", "Total"),values = c("#E69F00","#D55E00", "#0072B2", "#009E73", "#CC79A7"))
dev.off()
