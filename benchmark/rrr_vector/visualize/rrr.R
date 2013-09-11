require(tikzDevice)
source("../../basic_functions.R")

tex_file = "rrr.tex"

# Load experiment information
tc_config <- readConfig("../test_case.config",c("TC_ID","PATH","LATEX-NAME","URL"))
compile_config <- readConfig("../compile_options.config",c("COMPILE_ID","OPTIONS"))

open_tikz <- function( file_name  ){
    tikz(file_name, width = 5.5, height = 6, standAlone = F)
}

# Method which plots the space figure
plot_rrr_space <- function(data, max_y, title="", yaxis=T, xaxis=T, color){
  data <- data[order(data[['K']]), ]
  plot(c(), c(), ylim = c(0, max_y), xlim = c(min(data[['K']])-1, max(data[['K']])+1),
       yaxt = "n", ylab = "", xlab="", xaxt="n")

  if ( yaxis ){
      axis( 2, at = axTicks(2) , cex.axis=0.8)
      mtext("Space in (\\%) of original bitvector", side=2, line=2, las=0)
  }
  xticks <- c(8, 16, 32, 64, 126, 256)
  xlabel <- NA
  if ( xaxis ){
    xlabel <- c("8","","32","64","128","256")
    mtext("Block size K", side=1, line=2, las=0)
  }
  axis( 1, at = xticks, label=xlabel, cex.axis=0.8)
  hty =1 
  hcol = "lightgray"
  abline( v=c(8, 16, 32, 64, 128, 256), col=hcol )
  abline( h=axTicks(2), col=hcol )
  total_space <- 100*data[['rrr_size']]/data[['plain_size']][1]
  space_without_C <- 100*(data[['rrr_size']]-data[['bt_size']])/data[['plain_size']][1]
  space_only_samples <- 100*(data[['rrr_size']]-data[['bt_size']]-data[['btnr_size']])/data[['plain_size']][1]

  polygon( x_for_polygon( data[['K']] ), y_for_polygon( total_space ), border=NA, col=color[1])    
  polygon( x_for_polygon( data[['K']] ), y_for_polygon( space_without_C), border=NA, col=color[2])
  polygon( x_for_polygon( data[['K']] ), y_for_polygon( space_only_samples ), border=NA, col=color[3])
  lines( data[['K']], total_space )
  lines( data[['K']], space_without_C )
  lines( data[['K']], space_only_samples )
  draw_figure_heading(sprintf("bitvector = %s",title))
}

# Method which plots the query time figure
plot_rrr_query_times <- function( data, max_y=NA, title="", yaxis=T, xaxis=T){
  cat(title,"\n")
  data <- aggregate(data[c('access_time','rank_time','select_time')], by=c(data['K']), FUN=min)
  data[c('access_time','rank_time','select_time')] <- data[c('access_time','rank_time','select_time')]/1000.0
  data <- data[order(data[['K']]), ]

  max_runtime <- max( data[['access_time']], data[['rank_time']], data[['select_time']])
  if ( !is.na(max_y) ){
    max_runtime = max_y
  }

  plot(data[['K']], data[['access_time']], type="l", ylim=c(0, max_runtime ),
       xlab = "", ylab="", yaxt="n", cex.axis = 0.8, xaxt="n",
       col = terrain.colors(6)[1], lwd=1.5
       )
  box("plot", col="grey")    
  axis( 1, at = axTicks(1), labels=xaxis, mgp=c(2,0.5,0), tcl=-0.2, cex.axis=1, las=1 )
  if ( xaxis ){
    xlable <- "Block size K"
    mtext(xlable, side=1, line=2, las=0)
  }
  if ( yaxis ){
    axis( 2, at = axTicks(2),  mgp=c(1,0.3,0), tcl=-0.2, cex.axis=0.8 )
    mtext("Time per operation in ($\\mu s$)", side=2, line=2, las=0)
  }
  grid(lty=1)
  lines( data[['K']], data[['rank_time']], lty=1, col=terrain.colors(6)[3], lwd=1.5)    
  lines( data[['K']], data[['select_time']], lty=1, col=terrain.colors(6)[5], lwd=1.5)    

  draw_figure_heading(sprintf("bitvector = %s",title))
}

data <- data_frame_from_key_value_pairs( "../results/all.txt" )

tex_doc <- paste(readLines("rrr-header.tex"),collapse="\n")

for ( compile_id in compile_config[["COMPILE_ID"]] ){
# Handle query time
  fig_name <- paste("fig-rrr-time-",compile_id,".tex",sep="")
  open_tikz( fig_name )
  n <- nrow(tc_config)
  d <- subset(data, data[["COMPILE_ID"]]==compile_id)
  par(mfrow=c(n/2+1,2))
  multi_figure_style( n/2+1, 2 )  
  nr <- 0
  xlabnr <- 2*(n/2)-2+1
  for ( tc_id in tc_config[["TC_ID"]] ){
    dd <- subset(d, d[["TC_ID"]]==tc_id)
    plot_rrr_query_times(dd, 2.1, 
                         title=tc_config[tc_id,"LATEX-NAME"],
                            yaxis=(nr%%2==0), xaxis=(nr>=xlabnr) )
    if ( nr == 0 ){
        plot(NA, NA, xlim=c(0,1),ylim=c(0,1),ylab="", xlab="", bty="n", type="n", yaxt="n", xaxt="n")
        legend("topleft", legend=rev(c("access","rank","select")), box.lwd=0, lty=rev(c(1,1,1)), 
              title="Operation", col=rev(terrain.colors(6)[seq(1,5,2)]), bg="white", cex=1.5)
        nr <- nr+1
    }
    nr <-nr+1 
  }
  dev.off()
  tex_doc <- paste(tex_doc,"\\begin{figure}
               \\input{",fig_name,"}
               \\caption{Runtime of \\texttt{rrr\\_vector} dependent on block size K.
               Compile options:
               \\texttt{",gsub("_","\\\\_",compile_config[compile_id, "OPTIONS"]),"}.
               }
              \\end{figure}")
}

# Handle space
fig_name <- "fig-rrr-space.tex"
open_tikz( fig_name )
n <- nrow(tc_config)
d <- subset(data, data[["COMPILE_ID"]]==compile_config[1,1])
max_size <- 100*max(d[["rrr_size"]]/d[["plain_size"]])

multi_figure_style( n/2+1, 2 )  

nr <- 0
xlabnr <- 2*(n/2)-2+1
colors = c('gray80','gray50','gray30')
for ( tc_id in tc_config[["TC_ID"]] ){
  dd <- subset(d, d[["TC_ID"]]==tc_id)
  plot_rrr_space(dd, max_size, title=tc_config[tc_id,"LATEX-NAME"],
                      yaxis=(nr%%2==0), xaxis=(nr>=xlabnr), color=colors)
  if ( nr == 0 ){
    plot(NA, NA, xlim=c(0,1),ylim=c(0,1),ylab="", xlab="", bty="n", type="n", yaxt="n", xaxt="n")
    legend_text <- rev(c('pointers/samples','encoded block','block type'))
    legend( "top", legend=legend_text, fill=colors, bg="white", box.col="white", title="Space of",cex=1.5 )
    nr <- nr+1
  }
  nr <-nr+1 
}
dev.off()
tex_doc <- paste(tex_doc,"\\begin{figure}
                 \\input{",fig_name,"}
                 \\caption{Space of \\texttt{rrr\\_vector} dependent on block size K.}.
                 \\end{figure}")


tex_doc <- paste(tex_doc, readLines("rrr-footer.tex"),collapse="\n")
sink(tex_file)
cat(tex_doc)
sink(NULL)

## Method which plots the construction time figure
#plot_rrr_construction_times <- function( file_name, max_y=NA, draw_legend=T, title=""){
#  data <- data_frame_from_key_value_pairs( file_name )
#
#  data <- aggregate(data['construct_time'], by=c(data['K']), FUN=min)
#  data <- data[order(data[['K']]), ]
#
#  max_runtime <- max(data[['construct_time']])/1000
#  if ( !is.na(max_y) ){
#    max_runtime = max_y
#  }
#  plot(data[['K']], data[['construct_time']]/1000,
#       type="l", ylab="construction time in seconds",
#       xlab="Block size R", ylim=c(0, max_runtime ))
#  title( main=title, line=0 )
#}
#

