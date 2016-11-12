require(tikzDevice)
library(gdata)
library(xtable)
library(ggplot2)
library(reshape2)

source("../../basic_functions.R")
scale_data <- data_frame_from_key_value_pairs("../results/scale.log")
scale_data$constructs_space <- scale_data$constructs_space/1024/1024/1024
scale_data$`constructs_space(VMEM)` <- scale_data$`constructs_space(VMEM)`/1024/1024/1024
scale_data$construct_sort_duration <- scale_data$construct_sort_duration/1000
scale_data$construct_duration <- scale_data$construct_duration/1000
scale_data$buildvec_duration <- scale_data$buildvec_duration/1000
scale_data$subtree_constructor_duration <- scale_data$subtree_constructor_duration/1000
scale_data$compression_space <- scale_data$compression_space/1024/1024/1024
scale_data$`compression_space(VMEM)` <- scale_data$`compression_space(VMEM)`/1024/1024/1024
scale_data$constructor_call_duration <- scale_data$constructor_call_duration/1000
scale_data$comp_word_it <- scale_data$comp_word_it/1000
scale_data$frequency_encoding <- scale_data$frequency_encoding/1000
scale_data$dac_compression <- scale_data$dac_compression/1000

scale_2 <- subset(scale_data, select=c('constructor_call_duration', 'threads'))

constructor_call_duration <- melt(scale_2, id.var="threads")
ggplot(constructor_call_duration, aes(x = threads, y = value, fill = variable)) + 
  geom_bar(stat = "identity")

scale_3 <- subset(scale_data, select=c('construct_sort_duration', 'construct_duration', 'buildvec_duration', 'threads'))

stacked_construction <- melt(scale_3, id.var="threads")
ggplot(stacked_construction, aes(x = threads, y = value, fill = variable)) + 
  geom_bar(stat = "identity")

scale_4 <- subset(scale_data, select=c('comp_word_it', 'frequency_encoding', 'dac_compression', 'threads'))
compression_stacked <- melt(scale_4, id.var="threads")
ggplot(compression_stacked, aes(x = threads, y = value, fill = variable)) + 
  geom_bar(stat = "identity")

scale_5 <- subset(scale_data, select=c('compression_time', 'threads'))
compression <- melt(scale_5, id.var="threads")
ggplot(compression, aes(x = threads, y = value, fill = variable)) + 
  geom_bar(stat = "identity")

