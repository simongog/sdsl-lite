require(tikzDevice)
source("../../basic_functions.R")

tex_file = "gm.tex"

# Load experiment information
algo_config <- readConfig("../algorithms.config",c("ALGO_ID","LATEX-NAME","PCH","LTY","COL"))
pattern_config <- readConfig("../patterns.config",c("TC_ID","SP-LEN","GAP"))

open_tikz <- function( file_name  ){
    tikz(file_name, width = 5.5, height = 6, standAlone = F)
}

# Method which plots the query time figure
plot_gm_query_times <- function( data, title=""){
  cat("Graph: ",title,"\n")
  #data <- aggregate(data[c('total_time_mus')], by=c(data['PATT_SAMPLE']), FUN=min)
  #data[c('total_time_mus')] <- data[c('total_time_mus')]/1000.0
  #data <- data[order(data[['PATT_SAMPLE']]), ]

  runtime <- data[['total_time_mus']] # / data[['num_results']]
  max_runtime <- max(runtime[!is.infinite(runtime) & !is.nan(runtime)])
  max_sample <- length(unique(data[['PATT_SAMPLE']]))

  plot(NA, NA, type="l", xlim=c(0, max_sample ), ylim=c(0, max_runtime ),
       xlab = "", ylab="", yaxt="n", xaxt="n")

  ALGO_IDs <- unique(data$ALGO)
  for(algo in ALGO_IDs){
    d <- data[data$ALGO==algo,]
    lines(d[['PATT_SAMPLE']], d[['total_time_mus']] # / d[['num_results']]
          ,lwd=1, type="b", 
          pch=algo_config[algo,"PCH"], 
          lty=algo_config[algo,"LTY"],
          col=algo_config[algo,"COL"])
    box("plot", col="grey")    
    axis( 1, at = axTicks(1), labels=T, mgp=c(2,0.5,0), tcl=-0.2, cex.axis=1, las=1 )
  }

  mtext("Test case", side=1, line=2, las=0)
  axis( 2, at = axTicks(2),  mgp=c(1,0.3,0), tcl=-0.2, cex.axis=0.8 )
  mtext("Total query time in ($ms$)", side=2, line=2, las=0)
  #mtext("Average query time in ($\\mu s$)", side=2, line=2, las=0)

  grid(lty=1)  
  draw_figure_heading(sprintf("collection = %s",title))

  plot(NA, NA, xlim=c(0,1),ylim=c(0,1),ylab="", xlab="", bty="n", type="n", yaxt="n", xaxt="n")
  legend( "top", legend=algo_config[ALGO_IDs,"LATEX-NAME"], pch=algo_config[ALGO_IDs,"PCH"], col=algo_config[ALGO_IDs,"COL"],
            lty=algo_config[ALGO_IDs,"LTY"], bty="n", y.intersp=1.5, ncol=2, title="Algorithm", cex=1.2)
}

get_sa_range_for <- function(data, coll, pattern_id) {
  d <- data[data$PATT_SAMPLE==pattern_id,]
  numbers <- as.numeric(unlist(strsplit(paste(d[["info"]]),",")))
  numbers <- numbers[!is.na(numbers)]
  return(sum(numbers));
}

create_table_for <- function(data, coll, algo) {
  d <- data[data$ALGO==algo,]
  cat("Table: ",coll, " ", algo,"\n")
  # layout: rows=patterns, cols=values of interest

  PATT_IDs <- unique(d$PATT_SAMPLE)
  ROWS <- unique(pattern_config[PATT_IDs,"SP-LEN"])
  COLS <- unique(pattern_config[PATT_IDs,"GAP"])

  fig_name <- paste("fig-gm-time-",coll,"-",algo,"-table.tex",sep="")
  sink(fig_name)
  cat("\\begin{center}")
  cat("\\begin{tabular}{|r|", rep("r|", length(COLS)), "}\n", sep="")
  cat("\\hline\n")

  table = c("",ROWS)
            
  for(col in COLS){
    cfg_for_col <- pattern_config[pattern_config$GAP==col,]
    pids <- sapply(ROWS, function(row) { cfg_for_col[cfg_for_col$"SP-LEN"==row,][["TC_ID"]] })

    table = paste(table, 
             c(gsub(",","--",col),sapply(pids, function(patt_id)
                {
                    dd <- d[d$PATT_SAMPLE==patt_id,]
                    sa_range <- get_sa_range_for(data,coll,patt_id)
                    time <- dd[["mean_time_mus"]]
                    return(paste("$", sprintf("%.3f",time), "\\mu s$/$",sprintf("%.3f",time/sa_range),"\\mu s$", sep=""))
                })), sep=" & ")
  }

  cat(table, "", sep="\\\\\\hline\n")
  cat("\\end{tabular}\n")
  cat("\\end{center}")
  sink(NULL)
  
  return(paste("\\begin{figure}
                \\input{",fig_name,"}
                \\caption{Query times of algorithm \\texttt{",algo_config[algo,"LATEX-NAME"],"} on \\texttt{",coll,"}.}
                \\end{figure}", sep=""))
}

data <- data_frame_from_key_value_pairs( "../results/all.txt" )
tex_doc <- paste(readLines("gm-header.tex"),collapse="\n")

colls <- unique(data$COLL_ID)
n <- length(colls)
for(coll in colls){
  coll_name <- 
  
  # Graph for total query time
  fig_name <- paste("fig-gm-time-",coll,".tex",sep="")
  open_tikz( fig_name )

  multi_figure_style( 1, 1 )  

  d <- data[data$COLL_ID==coll,]
  plot_gm_query_times(d, title=coll)
  dev.off()
  tex_doc <- paste(tex_doc,"\\begin{figure}
               \\input{",fig_name,"}
               \\caption{Query time on \\texttt{",coll,"} depending on gap size.
               }
              \\end{figure}")
              
  # tables for mean query times and result numbers
  ALGO_IDs <- unique(d$ALGO)
  for(algo in ALGO_IDs){
    tex_doc <- paste(tex_doc,create_table_for(d, coll, algo))
  }

  tex_doc <- paste(tex_doc,"\\clearpage")
}

tex_doc <- paste(tex_doc, readLines("gm-footer.tex"),collapse="\n")
sink(tex_file)
cat(tex_doc)
sink(NULL)


