require(tikzDevice) 
source("../../basic_functions.R")


mypaste <- function(...){
    paste(...,sep="")
}

my_format2 <- function(...){
   format(..., nsmall=2,digits=2, big.mark=",")
}

process_data <- function(suf){

# Files
f_idxfilterconfig = mypaste("index-filter",suf,".config")
f_idxconfig       = mypaste("../index",suf,".config")
f_tcconfig        = mypaste("../test_case",suf,".config")
f_results         = mypaste("../results/all",suf,".txt")
f_sizes           = mypaste("../info/sizes",suf,".txt")
f_fig_runtime     = mypaste("fig-runtime",suf,".tex")
f_tbl_indexes     = mypaste("tbl-indexes",suf,".tex")
f_tbl_sizes       = mypaste("tbl-sizes",suf,".tex")
f_tbl_collections = mypaste("tbl-collections",suf,".tex")

# Load filter information
config <- readConfig(f_idxfilterconfig,c("IDX_ID","PCH","LTY","COL"))
# Load index and test case data
idx_config <- readConfig(f_idxconfig,c("IDX_ID","SDSL_TYPE","LATEX-NAME"))
tc_config <- readConfig(f_tcconfig,c("TC_ID","PATH","LATEX-NAME","URL"))

# Load data
raw <- data_frame_from_key_value_pairs(f_results)
raw <- raw[c("TC_ID","IDX_ID","time_per_query","query_len","TLE","doc_cnt","word_cnt")]
raw["time_per_query"] <- raw["time_per_query"]/1000.0

# Filter indexes
raw               <- raw[raw[["IDX_ID"]]%in%config[["IDX_ID"]],] 
raw[["IDX_ID"]]   <- factor(raw[["IDX_ID"]])

# Split by TC_ID
d <- split(raw,raw["TC_ID"])

tikz(f_fig_runtime, width = 5.5, height = 6, standAlone = F)

multi_figure_style( length(d)/2+1, 2 )  
count <- 0
nr <- 0
xlabnr <- (2*(length(d)/2)-1) 

for( tc_id in names(d) ){        

    plot( c(1e-9), c(1e-9), xlim=c(min(raw["query_len"]),max(raw["query_len"])),
                    ylim=c(0.01,1000),
                    ylab="", xlab="", xaxt="n", log="y" )
    box(col="gray")
    abline(h=c(seq(1,10)*0.1,seq(2,10)*1,seq(2,10)*10, seq(2,10)*100), col="lightgray")
    grid(lty="solid")
    if ( nr %% 2 == 0 ){
      ylable <- "Time per query (milliseconds)" 
      axis( 2, at = axTicks(2) )
      mtext(ylable, side=2, line=2, las=0)
    }
    axis( 1, at = axTicks(1), labels=(nr>=xlabnr)  )
    if ( nr >= xlabnr ){
      xlable <- "Pattern length"
      mtext(xlable, side=1, line=2, las=0)
    }
#   Split by IDX_ID
    dd <- split(d[[tc_id]],d[[tc_id]]["IDX_ID"])
    for( idx_id in names(dd) ){
        ddd <- dd[[idx_id]]
        ddd <- subset(ddd, ddd[["TLE"]]==0)
        lines(ddd[["query_len"]], ddd[["time_per_query"]],
              lwd=1, type="b", pch=config[idx_id, "PCH"], 
			  lty=config[idx_id, "LTY"],
			  col=config[idx_id, "COL"])
    }

    draw_figure_heading( sprintf("instance = \\textsc{%s}",tc_config[tc_id,"LATEX-NAME"]) )

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

sink(f_tbl_indexes)
cat(typeInfoTable(f_idxconfig, config[["IDX_ID"]], 1, 3, 2))
sink(NULL)

sink(f_tbl_sizes)
raw <- data_frame_from_key_value_pairs(f_sizes)
# Filter indexes
raw               <- raw[raw[["IDX_ID"]]%in%config[["IDX_ID"]],] 
raw[["IDX_ID"]]   <- factor(raw[["IDX_ID"]])

d <- split(raw,raw[["IDX_ID"]])

cat(paste("\\begin{tabular}{@{}l*{",length(names(d)) ,"}{rr@{\ }l}@{}}",sep=""))
cat("\\toprule\n")
cat(paste("Collection& & \\multicolumn{",length(names(d))*3-1,"}{c}{Index size in MiB (fraction of original collection)}\\\\\n",sep=""))
cat(paste("\\cmidrule{1-1}\\cmidrule{3-",length(names(d))*3+1,"}\n",sep=""))
for( idx_id in names(d) ){        
    cat(paste("&\\ &\\multicolumn{2}{c}{",idx_config[idx_id,"LATEX-NAME"],"}",sep=""))
}
cat("\\\\[2ex]\n")

dd <- split(raw,raw[["TC_ID"]])

for ( tc_id in names(dd) ){
    cat("{\\sc ",tc_config[tc_id,"LATEX-NAME"],"}")
    ddd <- split(dd[[tc_id]],dd[[tc_id]][["IDX_ID"]])
    for ( idx_id in names(d) ){
        row <- ddd[[idx_id]]
        cat(paste("&&",my_format2(row["size"]/1024**2),"&",sprintf("(%.2f)",row["size"]/row["text_size"])))
    }
    cat("\\\\\n")
}
cat("\\bottomrule\n")
cat("\\end{tabular}")

sink(NULL)

sink(f_tbl_collections)
raw <- data_frame_from_key_value_pairs(f_results)

cat(paste("\\begin{tabular}{@{}l*{5}{r}@{}}\n",sep=""))
cat("\\toprule\n")
elem = "Characters"
if ( suf=="_int" ) {
    elem = "Words"
}
cat(paste("Collection & ",elem," & Documents & ","Avg. doc. len.& gzip-compr.& xz-compr.\\\\[1ex]\n",sep=""))

d <- split(raw,raw["TC_ID"])
for ( tc_id in names(d) ){
    zfile = paste("../",tc_config[tc_id,"PATH"],".z.info",sep="")
    gzipcomp = NA
    xzcomp   = NA
    if ( file.exists(zfile) ){
        comp_info <- readConfig(zfile,c("NAME","RATIO","COMMAND"))
        gzipcomp = my_format2(comp_info["gzip","RATIO"])
        xzcomp   = my_format2(comp_info["xz","RATIO"])
    }
    doc_cnt_i  <- unique(d[[tc_id]][["doc_cnt"]])
    word_cnt_i <- unique(d[[tc_id]][["word_cnt"]])
    doc_cnt    <- format(doc_cnt_i, big.mark=",")
    word_cnt   <- format(word_cnt_i, big.mark=",")
    avg_doc_len <- my_format2(word_cnt_i/doc_cnt_i)
    cat("{\\sc ",tc_config[tc_id,"LATEX-NAME"],"}")
    cat("&",word_cnt,"&",doc_cnt,"&",avg_doc_len,"&",
        gzipcomp, "&", xzcomp,"\\\\\n")
}
cat("\\bottomrule\n")
cat("\\end{tabular}")
sink(NULL)

}


process_data("")
process_data("_int")
