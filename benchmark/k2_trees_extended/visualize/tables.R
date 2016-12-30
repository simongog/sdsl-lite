require(tikzDevice)
library(gdata)
library(xtable)
library(ggplot2)
library(reshape2)

preprocess_k2 <- function(df) {
  df$construction_space <-  pmax(df$constructs_space, df$compression_space)/1024/1024
  df$constructs_space <- NULL
  df$compression_space <- NULL
  df$construction_time <- df$constructs_time + df$compression_time
  df$constructs_time <- NULL
  df$compression_time <- NULL
  df <- split(df , f = df$TC_ID)
  df <- calculate_bpes(df)
  return(df)
}

calculate_bpes <- function(df) {
  df[["EU2005"]]$compressed_size <- df[["EU2005"]]$compressed_size*8 / 19235140
  df[["INDO"]]$compressed_size <- df[["INDO"]]$compressed_size*8 / 194109311
  df[["UK2002"]]$compressed_size <- df[["UK2002"]]$compressed_size*8 / 298113762
  df[["ARAB"]]$compressed_size <- df[["ARAB"]]$compressed_size*8 / 639999458
  df[["UK07-05"]]$compressed_size <- df[["UK07-05"]]$compressed_size*8 / 3738733648
  return(df)
}


source("../../basic_functions.R")
blnR <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/conni/line_chart.log")
parthR <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/conni/part_hybrid_rankv5.log")
hybridR <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/conni/hybrid_rankv5_no_short.log") 
webgraphR <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/conni/webgraph_tests_preview_with_llp.log")

bln <- subset(blnR, select=c('K2_ID','TC_ID', 'constructs_time', 'compression_time', 'constructs_space', 'compression_space', 'compressed_size', 'adj_time_comp', 'neighbors_time_comp', 'reverse_neighbors_time_comp'))
parth <- subset(parthR, select=c('K2_ID','TC_ID', 'constructs_time', 'compression_time', 'constructs_space', 'compression_space', 'compressed_size', 'adj_time_comp', 'neighbors_time_comp', 'reverse_neighbors_time_comp'))
#hybrid <- subset(hybridR, select=c('K2_TEX_NAME','TC_ID', 'constructs_time', 'compression_time', 'constructs_space', 'compression_space', 'compressed_size', 'adj_time_comp', 'neighbors_time_comp', 'reverse_neighbors_time_comp'))
#webgraph construction space (mapped?, add external construction?) adj_time
webgraph <- subset(webgraphR, select=c('config',  'TC_ID', 'construction_time', 'compressed_size', 'neighbors_time_comp', 'reverse_neighbors_time_comp'))

parth$neighbors_time_comp = parth$neighbors_time_comp/1000
parth$adj_time_comp = parth$adj_time_comp/1000
parth$reverse_neighbors_time_comp = parth$reverse_neighbors_time_comp/1000

#hybrid$neighbors_time_comp = hybrid$neighbors_time_comp/1000
#hybrid$adj_time_comp = hybrid$adj_time_comp/1000
#hybrid$reverse_neighbors_time_comp = hybrid$reverse_neighbors_time_comp/1000

webgraph$neighbors_time_comp = webgraph$neighbors_time_comp/1000
webgraph$reverse_neighbors_time_comp = webgraph$reverse_neighbors_time_comp/1000
webgraph$compressed_size = webgraph$compressed_size * 1024

parthp <- preprocess_k2(parth)
#hybridp <- preprocess_k2(hybrid)
blnp <- preprocess_k2(bln)
webgraphp <- split(webgraph , f = webgraph$TC_ID)
webgraphp <- calculate_bpes(webgraphp)

GRAPH = "EU2005"
PARTHID = "PartH4_5_2_8_S2"
BLNPID = "LADRA_5"
WGID1 = "eu2005w3m3"
WGID2 = "eu2005w70m300"
WGID3 = "eu2005w70m1000"

parthpEU <- parthp[[GRAPH]]
parthpEU <- parthpEU[parthpEU$K2_ID==PARTHID, ]

blnEU <- blnp[[GRAPH]]
blnEU <- blnEU[blnEU$K2_ID==BLNPID, ]

wg1 <- webgraphp[[GRAPH]]
wg1 <- wg1[wg1$config==WGID1, ]
wg2 <- webgraphp[[GRAPH]]
wg2 <- wg2[wg2$config==WGID2, ]
wg3 <- webgraphp[[GRAPH]]
wg3 <- wg3[wg3$config==WGID3, ]

wg <- rbind(wg1,wg2,wg3)
wg$construction_space <- 0
wg$adj_time_comp<- 0

k2 <- rbind(blnEU,parthpEU)
colnames(k2)[colnames(k2)=="K2_ID"] <- "config"

all <- rbind(wg,k2)
all$TC_ID <- NULL

all$construction_time <- format(all$construction_time, digits=0, nsmall=0)
all$compressed_size <- format(all$compressed_size, digits=2, nsmall=2)
all$adj_time_comp <- format(all$adj_time_comp, digits=2, nsmall=2)
all$neighbors_time_comp <- format(all$neighbors_time_comp, digits=2, nsmall=2)
all$reverse_neighbors_time_comp <- format(all$reverse_neighbors_time_comp, digits=2, nsmall=2)
all$construction_space <- format(all$construction_space, digits=0, nsmall=0)

colnames(all)[colnames(all)=="construction_time"] <- "construction time [s]"
colnames(all)[colnames(all)=="construction_space"] <- "construction space [Mb]"
colnames(all)[colnames(all)=="compressed_size"] <- "compressed size [bpe]"
colnames(all)[colnames(all)=="adj_time_comp"] <- "link check [$\\mu$s]"
colnames(all)[colnames(all)=="neighbors_time_comp"] <- "successors [$\\mu$s]"
colnames(all)[colnames(all)=="reverse_neighbors_time_comp"] <- "predecessors [$\\mu$s]"


xtable(t(all))

#sdslK2RF <- sdslk2R[grep("*_S", sdslk2R$TC_ID, invert=TRUE),]