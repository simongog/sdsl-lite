require(tikzDevice)
library(gdata)
library(xtable)
library(ggplot2)
library(reshape2)

lineplot <- function(name, direction, sdsl, ladra, webgraph){
  sdsl$TC_ID <- NULL
  sdsl <- melt(sdsl, id.var="compressed_size")
  sdsl$variable <- "SDSL"
  sdsl$value <- sdsl$value/1000   #sdsl ouptuts nsec
  
  webgraph$TC_ID <- NULL
  webgraph <- melt(webgraph, id.var="compressed_size")
  webgraph$variable <- "WEBGRAPH"
  webgraph$value <- webgraph$value/1000   #sdsl ouptuts nsec
  
  ladra$TC_ID <- NULL
  ladra <- melt(ladra, id.var="compressed_size")
  ladra$variable <- "BRIS"
  
  combined <- rbind(sdsl,ladra, webgraph)
  plotTitle <- paste(name, " - ", direction)
  #tikz(file = paste(direction,"-",name,"-space-time.tex",sep=""), width = 5, height = 3)
  plot <- ggplot(data=combined, aes(x=compressed_size, y=value, group=variable, shape=variable, colour=variable)) +
    geom_line() +
    labs(title = plotTitle, x = "$x$ = compressed size [bpe]", y = "$y$ = access time [$\\mu$sec]")+
    geom_point(size=2)+
    scale_colour_manual(name = '',
                        labels = c('BRIS','SDSL','BV(d+r)'),values = c("#E69F00","#009E73", "#0072B2"))+
    scale_shape_manual(name = '',
                        labels = c('BRIS', 
                                   'SDSL','BV(d+r)'),values = c(19,17,13)) #+
    #coord_cartesian(xlim = c(0, 5), ylim = c(0,4)) 
   #print(plot)
  #dev.off()
}


construction_time_plot <- function(name, direction, sdsl, ladra, webgraph, webgraph_direct){
  sdsl$TC_ID <- NULL
  sdsl <- melt(sdsl, id.var="compressed_size")
  sdsl$variable <- "SDSL"
  
  webgraph$TC_ID <- NULL
  webgraph <- melt(webgraph, id.var="compressed_size")
  webgraph$variable <- "WEBGRAPH"
  
  webgraph_direct$TC_ID <- NULL
  webgraph_direct <- melt(webgraph_direct, id.var="compressed_size")
  webgraph_direct$variable <- "WEBGRAPH-DIRECT"
  
  ladra$TC_ID <- NULL
  ladra <- melt(ladra, id.var="compressed_size")
  ladra$variable <- "BRIS"
  
  combined <- rbind(sdsl,ladra, webgraph, webgraph_direct)
  plotTitle <- paste(name, " - ", direction)
  tikz(file = paste(direction,"-",name,"-space-time.tex",sep=""), width = 5, height = 3)
  plot <- ggplot(data=combined, aes(x=compressed_size, y=value, group=variable, shape=variable, colour=variable)) +
    geom_line() +
    labs(title = plotTitle, x = "$x$ = compressed size [bpe]", y = "$y$ = construction time [sec]")+
    geom_point(size=4)+
    scale_colour_manual(name = '',
                        labels = c('BRIS','SDSL','BV(d+r)','BV(d)'),values = c("#E69F00","#009E73", "#0072B2","#D55E00"))+
    scale_shape_manual(name = '',
                       labels = c('BRIS', 
                                  'SDSL','BV(d+r)','BV(d)'),values = c(19,17,13,10))+
    scale_y_log10(breaks = c(100,200,500,1000,2000, 3000,5000,7000))
    #+
  #coord_cartesian(xlim = c(0, 5), ylim = c(0,4)) 
  print(plot)
  dev.off()
}


construction_time_access_time_plot <- function(name, direction, sdsl, ladra, webgraph, webgraph_direct){
  sdsl$TC_ID <- NULL
  sdsl <- melt(sdsl, id.var="neighbors_time_comp")
  sdsl$variable <- "SDSL"
  
  webgraph$TC_ID <- NULL
  webgraph <- melt(webgraph, id.var="neighbors_time_comp")
  webgraph$variable <- "WEBGRAPH"
  
  webgraph_direct$TC_ID <- NULL
  webgraph_direct <- melt(webgraph_direct, id.var="neighbors_time_comp")
  webgraph_direct$variable <- "WEBGRAPH-DIRECT"
  
  ladra$TC_ID <- NULL
  ladra <- melt(ladra, id.var="neighbors_time_comp")
  ladra$variable <- "BRIS"
  
  combined <- rbind(sdsl,ladra, webgraph, webgraph_direct)
  plotTitle <- paste(name, " - ", direction)
  #tikz(file = paste(direction,"-",name,"-space-time.tex",sep=""), width = 5, height = 3)
  plot <- ggplot(data=combined, aes(x=neighbors_time_comp, y=value, group=variable, shape=variable, colour=variable)) +
    geom_line() +
    labs(title = plotTitle, x = "$x$ = access time [bpe]", y = "$y$ = construction time [sec]")+
    geom_point(size=2)+
    scale_colour_manual(name = '',
                        labels = c('BRIS','SDSL','BV(d+r)','BV(d)'),values = c("#E69F00","#009E73", "#0072B2","#D55E00"))+
    scale_shape_manual(name = '',
                       labels = c('BRIS', 
                                  'SDSL','BV(d+r)','BV(d)'),values = c(19,17,13,10))+
    coord_cartesian(xlim = c(0, 5)) +
    scale_y_log10()
  
  print(plot)
  #dev.off()
}

source("../../basic_functions.R")
ladrak2 <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/conni/line_chart.log")
sdslk2 <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/conni/part_hybrid_rankv5.log")
sdslk2rank20 <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/conni/hybrid_parth/all.txt")
webgraph <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/conni/webgraph_tests_preview_with_llp.log")
test_data <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/test_data.data")

ladrak2R <- subset(ladrak2, select=c('TC_ID', 'compressed_size', 'reverse_neighbors_time_comp'))
sdslk2R <- subset(sdslk2, select=c('TC_ID','compressed_size', 'reverse_neighbors_time_comp'))
webgraphSub <- subset(webgraph, select=c('TC_ID', 'compressed_size', 'reverse_neighbors_time_comp'))

ladrak2R <- subset(ladrak2, select=c('TC_ID', 'compressed_size', 'constructs_time'))
sdslk2R <- subset(sdslk2, select=c('TC_ID','compressed_size', 'constructs_time'))
webgraphSub <- subset(webgraph, select=c('TC_ID', 'compressed_size', 'construction_time'))
webgraphDirect <- subset(webgraph, select=c('TC_ID', 'compressed_size_direct', 'construction_time_direct'))

ladrak2R <- subset(ladrak2, select=c('TC_ID', 'compressed_size', 'neighbors_time_comp'))
sdslk2R <- subset(sdslk2, select=c('TC_ID','compressed_size', 'neighbors_time_comp'))
webgraphSub <- subset(webgraph, select=c('TC_ID', 'compressed_size', 'neighbors_time_comp'))

ladrak2R <- subset(ladrak2, select=c('TC_ID', 'neighbors_time_comp', 'constructs_time'))
sdslk2R <- subset(sdslk2, select=c('TC_ID','neighbors_time_comp', 'constructs_time'))
sdslk2R$neighbors_time_comp = sdslk2R$neighbors_time_comp/1000
webgraphSub <- subset(webgraph, select=c('TC_ID', 'neighbors_time_comp', 'construction_time'))
webgraphSub$neighbors_time_comp = webgraphSub$neighbors_time_comp/1000
webgraphDirect <- subset(webgraph, select=c('TC_ID', 'neighbors_time_comp', 'construction_time_direct'))
webgraphDirect$neighbors_time_comp = webgraphDirect$neighbors_time_comp/1000


sdslK2RF <- sdslk2R[grep("*_S", sdslk2R$TC_ID, invert=TRUE),]
sEU <- sdslK2RF[sdslK2RF$TC_ID=="EU2005",]
sIndo <- sdslK2RF[sdslK2RF$TC_ID=="INDO",]
sUK02 <- sdslK2RF[sdslK2RF$TC_ID=="UK2002",]
sArab <- sdslK2RF[sdslK2RF$TC_ID=="ARAB",]
sUK07 <- sdslK2RF[sdslK2RF$TC_ID=="UK07-05",]
sEU$compressed_size = (sEU$compressed_size*8) / 19235140
sIndo$compressed_size = (sIndo$compressed_size*8) / 194109311
sUK02$compressed_size = (sUK02$compressed_size*8) / 298113762
sArab$compressed_size = (sArab$compressed_size*8) / 639999458
sUK07$compressed_size = (sUK07$compressed_size*8) / 3738733648


lEU <- ladrak2R[ladrak2R$TC_ID=="EU2005",]
lIndo <- ladrak2R[ladrak2R$TC_ID=="INDO",]
lUK02 <- ladrak2R[ladrak2R$TC_ID=="UK2002",]
lArab <- ladrak2R[ladrak2R$TC_ID=="ARAB",]
lUK07 <- ladrak2R[ladrak2R$TC_ID=="UK07-05",]
lEU$compressed_size = (lEU$compressed_size*8) / 19235140
lIndo$compressed_size = (lIndo$compressed_size*8) / 194109311
lUK02$compressed_size = (lUK02$compressed_size*8) / 298113762
lArab$compressed_size = (lArab$compressed_size*8) / 639999458
lUK07$compressed_size = (lUK07$compressed_size*8) / 3738733648

wEU <- webgraphSub[webgraphSub$TC_ID=="EU2005",]
wIndo <- webgraphSub[webgraphSub$TC_ID=="INDO",]
wUK02 <- webgraphSub[webgraphSub$TC_ID=="UK2002",]
wArab <- webgraphSub[webgraphSub$TC_ID=="ARAB",]
wUK07 <- webgraphSub[webgraphSub$TC_ID=="UK07-05",]
wEU$compressed_size = (wEU$compressed_size*8*1000) / 19235140
wIndo$compressed_size = (wIndo$compressed_size*8*1000) / 194109311
wUK02$compressed_size = (wUK02$compressed_size*8*1000) / 298113762
wArab$compressed_size = (wArab$compressed_size*8*1000) / 639999458
wUK07$compressed_size = (wUK07$compressed_size*8*1000) / 3738733648

names(webgraphDirect)[names(webgraphDirect) == 'compressed_size_direct'] <- 'compressed_size'
wdEU <- webgraphDirect[webgraphDirect$TC_ID=="EU2005",]
wdIndo <- webgraphDirect[webgraphDirect$TC_ID=="INDO",]
wdUK02 <- webgraphDirect[webgraphDirect$TC_ID=="UK2002",]
wdArab <- webgraphDirect[webgraphDirect$TC_ID=="ARAB",]
wdUK07 <- webgraphDirect[webgraphDirect$TC_ID=="UK07-05",]
wdEU$compressed_size = (wdEU$compressed_size*8*1000) / 19235140
wdIndo$compressed_size = (wdIndo$compressed_size*8*1000) / 194109311
wdUK02$compressed_size = (wdUK02$compressed_size*8*1000) / 298113762
wdArab$compressed_size = (wdArab$compressed_size*8*1000) / 639999458
wdUK07$compressed_size = (wdUK07$compressed_size*8*1000) / 3738733648

construction_time_access_time_plot("EU2005", "Construction", sEU, lEU, wEU, wdEU)
construction_time_access_time_plot("Indochina", "Construction", sIndo, lIndo, wIndo, wdIndo)
construction_time_access_time_plot("UK-2002", "Construction", sUK02, lUK02, wUK02, wdUK02)
construction_time_access_time_plot("Arabic", "Construction",sArab, lArab, wArab, wdArab)
construction_time_access_time_plot("UK-2007-05", "Construction" ,sUK07, lUK07, wUK07, wdUK07)

construction_time_plot("EU2005", "Construction", sEU, lEU, wEU, wdEU)
construction_time_plot("Indochina", "Construction", sIndo, lIndo, wIndo, wdIndo)
construction_time_plot("UK-2002", "Construction", sUK02, lUK02, wUK02, wdUK02)
construction_time_plot("Arabic", "Construction",sArab, lArab, wArab, wdArab)
construction_time_plot("UK-2007-05", "Construction" ,sUK07, lUK07, wUK07, wdUK07)

lineplot("EU2005", "Successors", sEU, lEU, wEU)
lineplot("Indochina", "Successors", sIndo, lIndo, wIndo)
lineplot("UK-2002", "Successors", sUK02, lUK02, wUK02)
lineplot("Arabic", "Successors",sArab, lArab, wArab)
lineplot("UK-2007-05", "Successors",sUK07, lUK07, wUK07)

lineplot("EU2005", "Predecessors", sEU, lEU, wEU)
lineplot("Indochina", "Predecessors", sIndo, lIndo, wIndo)
lineplot("UK-2002", "Predecessors", sUK02, lUK02, wUK02)
lineplot("Arabic", "Predecessors",sArab, lArab, wArab)
lineplot("UK-2007-05", "Predecessors",sUK07, lUK07, wUK07)
