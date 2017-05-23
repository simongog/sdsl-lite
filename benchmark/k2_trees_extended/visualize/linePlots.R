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
  ladra$variable <- "BLN"
  
  combined <- rbind(sdsl,ladra, webgraph)
  plotTitle <- paste(name, " - ", direction)
  #tikz(file = paste(direction,"-",name,"-space-time.tex",sep=""), width = 5, height = 3)
  plot <- ggplot(data=combined, aes(x=compressed_size, y=value, group=variable, shape=variable, colour=variable)) +
    geom_line() +
    labs(title = plotTitle, x = "$x$ = compressed size [bpe]", y = "$y$ = access time [$\\mu$sec]")+
    geom_point(size=2)+
    scale_colour_manual(name = '',
                        labels = c('BLN','SDSL','BV(d+r)'),values = c("#E69F00","#009E73", "#0072B2"))+
    scale_shape_manual(name = '',
                        labels = c('BLN', 
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
  ladra$variable <- "BLN"
  
  combined <- rbind(sdsl,ladra, webgraph, webgraph_direct)
  plotTitle <- paste(name, " - ", direction)
  tikz(file = paste(direction,"-",name,"-space-time.tex",sep=""), width = 5, height = 3)
  plot <- ggplot(data=combined, aes(x=compressed_size, y=value, group=variable, shape=variable, colour=variable)) +
    geom_line() +
    labs(title = plotTitle, x = "$x$ = compressed size [bpe]", y = "$y$ = construction time [sec]")+
    geom_point(size=4)+
    scale_colour_manual(name = '',
                        labels = c('BLN','SDSL','BV(d+r)','BV(d)'),values = c("#E69F00","#009E73", "#0072B2","#D55E00"))+
    scale_shape_manual(name = '',
                       labels = c('BLN', 
                                  'SDSL','BV(d+r)','BV(d)'),values = c(19,17,13,10))+
    scale_y_log10(breaks = c(100,200,500,1000,2000, 3000,5000,7000))
    #+
  #coord_cartesian(xlim = c(0, 5), ylim = c(0,4)) 
  print(plot)
  dev.off()
}

clearTC_ID_melt_name <- function(df, id_var, name){
  df$TC_ID <- NULL
  df <- melt(df, id.var=id_var)
  df$variable <- name
  df
}

construction_time_space_plot <- function(name, ladra, count, countp, szord, szordp, pzord, pzordp){
  ladra <- clearTC_ID_melt_name(ladra,'constructs_space', "BLN")
  count <- clearTC_ID_melt_name(count, "constructs_space", "COUNT")
  countp <- clearTC_ID_melt_name(countp, "constructs_space", "COUNTP")
  szord <- clearTC_ID_melt_name(szord, "constructs_space", "SZORD")
  szordp <- clearTC_ID_melt_name(szordp, "constructs_space", "SZORDP")
  pzord <- clearTC_ID_melt_name(pzord, "constructs_space", "PZORD")
  pzordp <- clearTC_ID_melt_name(pzordp, "constructs_space", "PZORDP")
            
  combined <- rbind(ladra, count, countp, szord, szordp, pzord, pzordp)
  
  plotTitle <- paste(name)
  tikz(file = paste(name,"-time-mem.tex",sep=""), width = 5, height = 3)
  plot <- ggplot(data=combined, aes(x=constructs_space, y=value, group=variable, shape=variable, colour=variable)) +
    geom_line() +
    labs(title = plotTitle, x = "$x$ = construction space [Mbyte]", y = "$y$ = construction time [sec]")+
    geom_point(size=4) +
    scale_colour_manual(name = '',
                        labels = c('BLNP','COUNT','COUNTP','PZORD','PZORDP','SZORD','SZORDP'),values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
    scale_shape_manual(name = '',
                       labels = c('BLNP','COUNT','COUNTP','PZORD','PZORDP','SZORD','SZORDP'),values = c(16,17,18,19,20,7,11))
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
  ladra$variable <- "BLN"
  
  combined <- rbind(sdsl,ladra, webgraph, webgraph_direct)
  plotTitle <- paste(name, " - ", direction)
  #tikz(file = paste(direction,"-",name,"-space-time.tex",sep=""), width = 5, height = 3)
  plot <- ggplot(data=combined, aes(x=neighbors_time_comp, y=value, group=variable, shape=variable, colour=variable)) +
    geom_line() +
    labs(title = plotTitle, x = "$x$ = access time [bpe]", y = "$y$ = construction time [sec]")+
    geom_point(size=2)+
    scale_colour_manual(name = '',
                        labels = c('BLN','SDSL','BV(d+r)','BV(d)'),values = c("#E69F00","#009E73", "#0072B2","#D55E00"))+
    scale_shape_manual(name = '',
                       labels = c('BLN', 
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

tmp <- sEU
tmp$TC_ID <- NULL
tmp <- melt(tmp, id.var="compressed_size")
d <- tmp[ order(tmp$compressed_size, decreasing=TRUE), ]
result <- d[1,]
for(i in seq_len(nrow(d))[-1] ) {
  if( d$value[i] > result$value[nrow(result)] ) {
    result <- rbind(result, d[i,])  # inefficient
  } 
}
points(result, cex=3, pch=15)

require(rPref)


#construction time + space
ladrak2 <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/conni/line_chart.log")
counting_sort <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/conni/counting_sort.log")
counting_sort_parth <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/conni/counting_sort_partitioned_hybrid.log")
seq_zord_parth <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/conni/sequential_z_order_sort_part_tree_seq_access.log")
seq_zord <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/conni/sequential_z_order_sort_tree_seq_access.log")
hybrid_rankv5 <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/conni/hybrid_rankv5_no_short.log") 
part_hybrid_rankv5 <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/conni/part_hybrid_rankv5.log")


#this is all without dac compression and data is with parallel dac compression!
ladra <- subset(ladrak2, select=c('TC_ID', 'constructs_time', 'constructs_space'))
csort <- subset(counting_sort, select=c('TC_ID', 'constructs_time', 'constructs_space'))
csortpart <- subset(counting_sort_parth, select=c('TC_ID', 'constructs_time', 'constructs_space'))
s_zord <- subset(seq_zord, select=c('TC_ID', 'constructs_time', 'constructs_space'))
s_zord_part <- subset(seq_zord_parth, select=c('TC_ID', 'constructs_time', 'constructs_space'))
p_zord <- subset(hybrid_rankv5, select=c('TC_ID', 'constructs_time', 'constructs_space'))
p_zord_part <- subset(part_hybrid_rankv5, select=c('TC_ID', 'constructs_time', 'constructs_space'))


ladra <- preprocess_construction(ladra)
csort <- preprocess_construction(csort)
csortpart <- preprocess_construction(csortpart)
s_zord <- preprocess_construction(s_zord)
s_zord_part <- preprocess_construction(s_zord_part)
p_zord <- preprocess_construction(p_zord)
p_zord_part <- preprocess_construction(p_zord_part)

name = "UK07-05"
construction_time_space_plot(name, ladra[[name]], csort[[name]], csortpart[[name]], s_zord[[name]], s_zord_part[[name]], p_zord[[name]], p_zord_part[[name]])




#csort <-subset(counting_sort, select=c('TC_ID','K2_TEX_NAME','TC_TEX_NAME', 'constructs_time', 'constructs_space', 'compression_time', 'compressed_size', 'adj_time_comp', 'neighbors_time_comp', 'reverse_neighbors_time_comp'))
#hybrid <-subset(hybrid_rankv5, select=c('TC_ID','K2_TEX_NAME','TC_TEX_NAME', 'constructs_time', 'constructs_space', 'compression_time', 'compressed_size', 'adj_time_comp', 'neighbors_time_comp', 'reverse_neighbors_time_comp'))
#parth <-subset(part_hybrid_rankv5, select=c('TC_ID','K2_TEX_NAME','TC_TEX_NAME', 'constructs_time', 'constructs_space', 'compression_time', 'compressed_size', 'adj_time_comp', 'neighbors_time_comp', 'reverse_neighbors_time_comp'))

csort$constructs_time <- csort$constructs_time + csort$compression_time
csort$compression_time <- NULL
hybrid$constructs_time <- hybrid$constructs_time + hybrid$compression_time
hybrid$compression_time <- NULL
parth$constructs_time <- parth$constructs_time + parth$compression_time
parth$compression_time <- NULL

preprocess_construction <- function(df) {
  df$constructs_space <-  df$constructs_space/1024/1024
  split(df , f = df$TC_ID)
}

csortEU <- csort[csort$TC_ID=="EU2005",]
hybridEu <- hybrid[hybrid$TC_ID=="EU2005",]
partEu <- parth[parth$TC_ID=="EU2005",]

csort4628 <- subset(csortEU[csortEU$K2_TEX_NAME=="Hybrid(4628)S2",], select=c('constructs_time', 'constructs_space', 'compressed_size', 'adj_time_comp', 'neighbors_time_comp', 'reverse_neighbors_time_comp'))
hybrid4628 <- subset(hybridEu[hybridEu$K2_TEX_NAME=="Hybrid(4628)S2",], select=c('constructs_time', 'constructs_space', 'compressed_size', 'adj_time_comp', 'neighbors_time_comp', 'reverse_neighbors_time_comp'))
part4628 <- subset(partEu[partEu$K2_TEX_NAME=="PartHybrid(4628)S2",], select=c('constructs_time', 'constructs_space', 'compressed_size', 'adj_time_comp', 'neighbors_time_comp', 'reverse_neighbors_time_comp'))

csort4628$construction <- "Counting Sort"
hybrid4628$construction <- "ParZord Hybrid"
part4628$construction <- "ParZord PartH"

combined <- rbind(csort4628, hybrid4628, part4628)
transposed <- t(combined)

combEU <- combined[combined$TC_ID=="EU2005",]

p <- low(constructs_time)*low(constructs_space)*low(compressed_size)*low(neighbors_time_comp, df = combEU)
optimal <- peval(p)

show_front <- function(wEu, pref) {
  plot(wEu$compressed_size, wEu$reverse_neighbors_time_comp)
  sky <- psel(wEu, pref)
  plot_front(wEu, pref, col = rgb(0, 0, 1))
  points(sky$compressed_size, sky$reverse_neighbors_time_comp, lwd = 3)
}


pieEu4628 <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/pie-charts/eu4628.log")
pieEu4628Comp <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/pie-charts/eu4628comp.log")

pieEuK4 <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/pie-charts/euk4.log")
pieEuK4Comp <- data_frame_from_key_value_pairs("/home/d056848/Dev/master-thesis/benchmarks/pie-charts/euk4comp.log")

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
library(scales)

pieEuK4Comp$size = pieEuK4Comp$size/sum(pieEuK4Comp$size)*100
pieEuK4$size = pieEuK4$size/sum(pieEuK4$size)*100



size <- c(92.5, 2.96, 1.51, 3.03)
name <- c("last level (11)", "level 10", "level 9", "rest")
eu4628 <- data.frame(size,name)

size <- c(56.4, 18.9, 9.7, 4.9, 4.9, 2.5, 2.7)
name <- c("last level (dac)", "dictionary", "level 10",  "level 9", "rank helper", "level 8", "rest")
eu4628Comp <- data.frame(size,name)

size <- c(68.6, 18.8, 6.3, 4.8, 1.5)
name <- c("last level (9)", "level 8", "rank helper",  "level 7", "rest")
euK4 <- data.frame(size,name)

size <- c(50.9, 28.9, 9.6, 7.3, 3.3)
name <- c("last level (dac)", "level 8", "rank helper",  "level 7", "rest")
euK4Comp <- data.frame(size,name)

tikz(file = "pie-k4comp.tex", width = 5, height = 3)
plot <- ggplot(euK4Comp, aes(x="", y=size, fill=name))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme+
  theme(axis.text.x=element_blank()) +
  geom_text(aes(x=c(1,1,1,1.2,1.4), y = size / 2 + c(0, cumsum(size)[-length(size)]),
    #y = size/3 + c(0, cumsum(size)[-length(size)]), 
                label = paste(size),"\\%",sep=""), size=5)
print(plot)
dev.off()

#1,1,1,1.1,1.2,1.3,1.4
