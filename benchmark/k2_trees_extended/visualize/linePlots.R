require(tikzDevice)
library(gdata)
library(xtable)
library(ggplot2)
library(reshape2)

lineplot <- function(sdsl, ladra){
  sdsl$TC_ID <- NULL
  sdsl <- melt(sdsl, id.var="compressed_size")
  sdsl$variable <- "SDSL"
  sdsl$value <- sdsl$value/1000   #sdsl ouptuts nsec
  
  ladra$TC_ID <- NULL
  ladra <- melt(ladra, id.var="compressed_size")
  ladra$variable <- "BRIS"
  
  combined <- rbind(sdsl,ladra)
  ggplot(data=combined, aes(x=compressed_size, y=value, group=variable, shape=variable, colour=variable)) +
    geom_line() +
    xlab("$x$ = compressed size [bpe]")+
    ylab("$y$ = access time [msec]")+
    geom_point(size=4)+
    scale_colour_manual(name = '',
                        labels = c('BRIS','SDSL'),values = c("#E69F00","#009E73"))+
    scale_shape_manual(name = '',
                        labels = c('BRIS', 
                                   'SDSL'),values = c(19,17))
  
}

source("../../basic_functions.R")
ladrak2 <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/conni/table3.log")
sdslk2 <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/conni/hybrid_parth/parth_all.txt")
test_data <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/test_data.data")

ladrak2R <- subset(ladrak2, select=c('TC_ID','TC_TEX_NAME', 'compressed_size', 'neighbors_time_comp', 'reverse_neighbors_time_comp',  'adj_time_comp')) #, 'K2_ID','K2_TEX_NAME','K2_TEX_NAME',))
sdslk2R <- subset(sdslk2, select=c('TC_ID','TC_TEX_NAME', 'compressed_size', 'neighbors_time_comp', 'reverse_neighbors_time_comp',  'adj_time_comp'))

ladrak2R <- subset(ladrak2, select=c('TC_ID', 'compressed_size', 'neighbors_time_comp'))
sdslk2R <- subset(sdslk2, select=c('TC_ID','compressed_size', 'neighbors_time_comp'))

sdslK2RF <- sdslk2R[grep("*_S", sdslk2R$TC_ID, invert=TRUE),]
sEU <- sdslK2RF[sdslK2RF$TC_ID=="EU2005",]
sEU$compressed_size = (sEU$compressed_size*8) / 19235140
sIndo <- sdslK2RF[sdslK2RF$TC_ID=="INDO",]
sIndo$compressed_size = (sIndo$compressed_size*8) / 194109311
sUK02 <- sdslK2RF[sdslK2RF$TC_ID=="UK2002",]
sUK02$compressed_size = (sUK02$compressed_size*8) / 298113762
sArab <- sdslK2RF[sdslK2RF$TC_ID=="ARAB",]
sArab$compressed_size = (sArab$compressed_size*8) / 639999458
sUK07 <- sdslK2RF[sdslK2RF$TC_ID=="UK07-05",]
sUK07$compressed_size = (sUK07$compressed_size*8) / 3738733648

lEU <- ladrak2R[ladrak2R$TC_ID=="EU2005",]
lEU$compressed_size = (lEU$compressed_size*8) / 19235140
lIndo <- ladrak2R[ladrak2R$TC_ID=="INDO",]
lIndo$compressed_size = (lIndo$compressed_size*8) / 194109311
lUK02 <- ladrak2R[ladrak2R$TC_ID=="UK2002",]
lUK02$compressed_size = (lUK02$compressed_size*8) / 298113762
lArab <- ladrak2R[ladrak2R$TC_ID=="ARAB",]
lArab$compressed_size = (lArab$compressed_size*8) / 639999458
lUK07 <- ladrak2R[ladrak2R$TC_ID=="UK07-05",]
lUK07$compressed_size = (lUK07$compressed_size*8) / 3738733648

lineplot(sEU, lEU)
lineplot(sIndo, lIndo)
lineplot(sUK02, lUK02)
lineplot(sArab, lArab)
lineplot(sUK07, lUK07)
