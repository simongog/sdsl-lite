require(tikzDevice)
source("../../basic_functions.R")

# Load filter information
config <- readConfig("index-filter.config",c("IDX_ID","PCH","LTY","COL"))
idx_config <- readConfig("../index.config",c("IDX_ID","SDSL_TYPE","LATEX-NAME"))
tc_config <- readConfig("../test_case.config",c("TC_ID","PATH","LATEX-NAME","URL"))

# Load data
raw <- data_frame_from_key_value_pairs( "../results/all.txt" )
# Filter indexes
raw               <- raw[raw[["IDX_ID"]]%in%config[["IDX_ID"]],] 
raw[["IDX_ID"]] <- factor(raw[["IDX_ID"]])
# Normalize data
raw[["Time"]] <- 1000000*raw[["Extract_time_in_sec"]]/raw[["Total_num_chars_extracted"]]
raw[["Space"]]    <- 100*raw[["Index_size_in_bytes"]]/raw[["text_size"]]

raw <- raw[c("TC_ID", "Space", "Time","IDX_ID","S_SA","S_ISA")]
raw <- raw[order(raw[["TC_ID"]]),]

data <- split(raw, raw[["TC_ID"]])

tikz("fig-extract.tex", width = 5.5, height = 6, standAlone = F)

multi_figure_style( length(data)/2+1, 2 )  

max_space <- 100*2.5 #1.3 # 2.5 #1.8 
max_time  <- max(raw[["Time"]]) # in microseconds
count <- 0
nr <- 0
xlabnr <- (2*(length(data)/2)-1) 
for( tc_id in names(data) ){
  d <- data[[tc_id]]

  plot(c(),c(),xlim=c(0, max_space), ylim=c(0, max_time), xlab="", axes=F, xaxt="n", yaxt="n", ylab="" )
  box(col="gray")
  grid(lty="solid")
  if ( nr %% 2 == 0 ){
    ylable <- "Time per character ($\\mu s$)" 
    axis( 2, at = axTicks(2) )
    mtext(ylable, side=2, line=2, las=0)
  }
  axis( 1, at = axTicks(1), labels=(nr>=xlabnr) )
  if ( nr >= xlabnr ){
    xlable <- "Index size in (\\%)"
    mtext(xlable, side=1, line=2, las=0)
  }
  dd <- split(d, d[["IDX_ID"]])
  for( idx_id in names(dd) ){
    ddd <- dd[[idx_id]]
    lines(ddd[["Space"]], ddd[["Time"]], type="b", lwd=1, pch=config[idx_id, "PCH"], 
                                                          lty=config[idx_id, "LTY"],
                                                          col=config[idx_id, "COL"]
         )
  }
  draw_figure_heading( sprintf("instance = \\textsc{%s}",tc_config[tc_id,"LATEX-NAME"]) )
  comp_info <- readConfig(paste("../",tc_config[tc_id,"PATH"],".z.info",sep=""),c("NAME","RATIO","COMMAND"))
  abline(v=comp_info["xz","RATIO"],lty=1);
  abline(v=comp_info["gzip","RATIO"],lty=4);

  nr <- nr+1
  if ( nr == 1 ){ # plot legend
    plot(NA, NA, xlim=c(0,1),ylim=c(0,1),ylab="", xlab="", bty="n", type="n", yaxt="n", xaxt="n")
    idx_ids <- as.character(unique(raw[["IDX_ID"]]))
    legend( "top", legend=idx_config[idx_ids,"LATEX-NAME"], pch=config[idx_ids,"PCH"], col=config[idx_ids,"COL"],
            lty=config[idx_ids,"LTY"], bty="n", y.intersp=1.5, ncol=2, title="Index", cex=1.2)
    legend( "bottom", legend=c("\\textsc{xz}~(\\texttt{-9})","\\textsc{gzip}~(\\texttt{-9})"), 
            lty=c(1,4), bty="n", y.intersp=1.5, ncol=2, title="Compression baseline", cex=1.2)
    nr <- nr+1
  }
}
dev.off()

sink("tbl-extract.tex")
cat(typeInfoTable("../index.config", config[["IDX_ID"]], 1, 3, 2))
sink(NULL)
