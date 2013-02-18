library(xtable) # if not installed call install.packages("xtable")

source("basic_functions.R")

data_sdsl<- data_frame_from_key_value_pairs( "../results/count.txt" )
data_sdsl <- data_sdsl[data_sdsl[["program"]]%in%c("FM_HUFF","FM_HUFF_RRR63","FM_RLMN","FM_HUFF_RRR15","FM_HUFF_RRR127","FM_HUFF_RRR255"),] 
data_sdsl[["program"]] <- factor(data_sdsl[["program"]])

data_sdsl[["Time"]] <- data_sdsl[["Count_time_in_milli_sec"]]*1000/(data_sdsl[["pat_cnt"]]*data_sdsl[["pat_length"]])
data_sdsl[["Time"]] <- round(data_sdsl[["Time"]], 3)

data_sdsl[["Space"]] <- data_sdsl[["Index_size_in_bytes"]]/data_sdsl[["text_size"]]
data_sdsl[["Space"]] <- round(data_sdsl[["Space"]],2)

d2 <- data_sdsl[c("program","test_case", "Space", "Time","hugepages","sse","popcount_tl")]
d2 <- d2[order(d2[["test_case"]]),]

d_NOOPT <- subset(d2, d2["hugepages"]==0 & d2["sse"]==0 & d2["popcount_tl"]==1 )
d_NOSSE <- subset(d2, d2["hugepages"]==0 & d2["sse"]==0 & d2["popcount_tl"]==0 )
d_SSE <- subset(d2, d2["hugepages"]==0 & d2["sse"]==1 & d2["popcount_tl"]==0 )
d_HP <- subset(d2, d2["hugepages"]==1 & d2["sse"]==1 & d2["popcount_tl"]==0 )

form_table <- function(d, order=NA){
	d <- aggregate(d[c('Time','Space')], 
			       by=list(program=d[['program']],test_case=d[['test_case']],
					       hugepages=d[['hugepages']],sse=d[['sse']]),
				   FUN=mean,na.rm=TRUE)
	dByProgram <- split(d, d[["program"]])
	table <- data.frame(dByProgram[[1]]["test_case"])
	names(table) <- c("Text")
	prog_name <- names(dByProgram)
	if( !is.na(order) ){
		prog_name <- order
	}
	for( prog in prog_name ){
		sel <- dByProgram[[prog]]
		table <- cbind(table, sel["Time"])
		table <- cbind(table, sel["Space"])
	}
	list("table" = table, "names" = prog_name)
}


make_latex_header <- function(names){
	x <- paste("&\\multicolumn{2}{c|}{", names,"}")
	x <- paste(x, collapse=" ")
	y <- paste("\\hline",x, "\\\\\\hline\n")	
	gsub("_","\\\\_",y)
}

print_latex <- function( table, names ){
	ali <- c("|l|","|l|", rep("c|", ncol(table)-1)) 
	dig <- c(0, 0, rep(c(3,2),(ncol(table)-1)/2 ))	
	print( xtable( table, align=ali, digits=dig ), 
		   type="latex",  hline.after = c(0,0,seq(1, nrow(table))),
		   floating = F, # don't use table environment
		   add.to.row=list(pos=list(-1), 
				           command=(make_latex_header(names))),
		   sanitize.rownames.function = identity, 
		   sanitize.text.function = identity,
		   include.rownames = FALSE
		 )
}

sink("tbl-count-NOOPT.tex")
x3 <- form_table(d_NOOPT,c("FM_HUFF","FM_HUFF_RRR15","FM_RLMN","FM_HUFF_RRR63"))
print_latex(x3[["table"]], x3[["names"]])	
sink(NULL)

sink("tbl-count-NOSSE.tex")
x4 <- form_table(d_NOSSE,c("FM_HUFF","FM_HUFF_RRR15","FM_RLMN","FM_HUFF_RRR63"))
print_latex(x4[["table"]], x4[["names"]])	
sink(NULL)

sink("tbl-count-SSE.tex")
x5 <- form_table(d_SSE,c("FM_HUFF","FM_HUFF_RRR15","FM_RLMN","FM_HUFF_RRR63"))
print_latex(x5[["table"]], x5[["names"]])	
sink(NULL)

sink("tbl-count-HP.tex")
if ( nrow(d_HP) > 0 ){
	x6 <- form_table(d_HP,c("FM_HUFF","FM_HUFF_RRR15","FM_RLMN","FM_HUFF_RRR63"))
	print_latex(x6[["table"]], x6[["names"]])	
}else{
	cat("\\begin{center}1 GB pages were not supported\\end{center}")
}
sink(NULL)
