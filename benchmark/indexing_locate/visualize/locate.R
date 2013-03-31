if ( !exists( "tikzDeviceLoaded" ) ){  
	require(tikzDevice) #if not installed call install.packages("tikzDevice", repos="http://R-Forge.R-project.org")
	tikzDeviceLoaded = T
}
source("../../basic_functions.R")

# Load filter information
config <- readConfig("index-filter.config",c("IDX_ID","PCH","LTY","COL"))
idx_config <- readConfig("../index.config",c("IDX_ID","SDSL_TYPE","LATEX-NAME"))
tc_config <- readConfig("../test_case.config",c("TC_ID","PATH","LATEX-NAME"))

# Load data
raw <- data_frame_from_key_value_pairs( "../results/all.txt" )
# Filter indexes
raw               <- raw[raw[["IDX_ID"]]%in%config[["IDX_ID"]],] 
raw[["IDX_ID"]] <- factor(raw[["IDX_ID"]])
# Normalize data
raw[["Time"]] <- 1000000*raw[["Locate_time_in_secs"]]/raw[["Total_Num_occs_found"]]
raw[["Space"]]    <- 100*raw[["Index_size_in_bytes"]]/raw[["text_size"]]

raw <- raw[c("TC_ID", "Space", "Time","IDX_ID","S_SA","S_ISA")]
raw <- raw[order(raw[["TC_ID"]]),]

data <- split(raw, raw[["TC_ID"]])

tikz("locate.tex", width = 5.5, height = 6, standAlone = T)

par(mfrow=c((length(data)+2)/2, 2))

par(oma=c(2.5,2.7,0,0))
par(mar=c(1,1,1.5,0.5))
par(mgp=c(2,0.5,0))	

par(yaxs="i") # don't add +- 4% to the yaxis
par(xaxs="i") # don't add +- 4% to the xaxis

max_space <- 100*2.5 #1.3 # 2.5 #1.8 
max_time  <- 25    # in microseconds
count <- 0
nr <- 0
for( tc_id in names(data) ){
  d <- data[[tc_id]]

  plot(c(),c(),xlim=c(0, max_space), ylim=c(0, max_time), xlab="", axes=F, xaxt="n", yaxt="n", ylab="" )
  box(col="gray")
  grid(lty="solid")
  if ( nr %% 2 == 0 ){
    ylable <- "Time per occurrence ($\\mu s$)" 
    axis( 2, at = axTicks(2),  mgp=c(2,0.5,0), tcl=-0.2, cex.axis=1, las=1 )
    mtext(ylable, side=2, line=2, las=0)
  }
  axis( 1, at = axTicks(1), labels=(nr>3), mgp=c(2,0.5,0), tcl=-0.2, cex.axis=1, las=1 )
  if ( nr >= (2*(length(data)/2)-1) ){
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
  rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[4], ytop=par("usr")[4]*1.10 ,xpd=NA,
       col="grey80", border="grey80" )
  text(labels=sprintf("instance = \\textsc{%s}",tc_config[tc_id,"LATEX-NAME"]),y=par("usr")[4]*1.03,adj=c(0.5, 0),x=(par("usr")[1]+par("usr")[2])/2,xpd=NA,cex=1.4)
  nr <- nr+1
  if ( nr == 1 ){ # plot legend
    plot(NA, NA, xlim=c(0,1),ylim=c(0,1),ylab="", xlab="", bty="n", type="n", yaxt="n", xaxt="n")
	idx_ids <- as.character(unique(raw[["IDX_ID"]]))
    legend( "top", legend=idx_config[idx_ids,"LATEX-NAME"], pch=config[idx_ids,"PCH"], col=config[idx_ids,"COL"],
		    lty=config[idx_ids,"LTY"], bty="n", y.intersp=1.5, ncol=2, title="Index", cex=1.2)
    nr <- nr+1
  }
}
dev.off()
