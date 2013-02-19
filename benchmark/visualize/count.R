library(xtable) # if not installed call install.packages("xtable")

source("basic_functions.R")

#load index information
index_info <- read.csv("../index.config", sep=";",header=F, comment.char="#")
colnames(index_info)<-c("id","sdsl_type","latex_string")
id2latex <- data.frame( t(as.character(index_info[[3]])) )
names(id2latex) <- as.character(index_info[[1]]) 

#load report information
report_info <- read.csv("count.config", sep=";",header=F, comment.char="#")
index_ids <- as.character(unlist(report_info))

# tranforms a vector of index ids to a vector which contains the
# corresponding latex names.
# Note: each id should only appear once in the input vector
index_id2latex_name <- function(ids){
	as.character(unlist(id2latex[ids]))
}

data_sdsl<- data_frame_from_key_value_pairs( "../results/count.txt" )
data_sdsl <- data_sdsl[data_sdsl[["program"]]%in%index_ids,] 
data_sdsl[["program"]] <- factor(data_sdsl[["program"]])

data_sdsl[["Time"]] <- data_sdsl[["Count_time_in_milli_sec"]]*1000/(data_sdsl[["pat_cnt"]]*data_sdsl[["pat_length"]])
data_sdsl[["Time"]] <- round(data_sdsl[["Time"]], 3)

data_sdsl[["Space"]] <- data_sdsl[["Index_size_in_bytes"]]/data_sdsl[["text_size"]]
data_sdsl[["Space"]] <- round(data_sdsl[["Space"]],2)

d2 <- data_sdsl[c("program","test_case", "Space", "Time","hugepages","sse","popcount_tl")]
d2 <- d2[order(d2[["test_case"]]),]

data <- list()
data[["NOOPT"]] <- subset(d2, d2["hugepages"]==0 & d2["sse"]==0 & d2["popcount_tl"]==1 )
data[["NOSSE"]] <- subset(d2, d2["hugepages"]==0 & d2["sse"]==0 & d2["popcount_tl"]==0 )
data[["SSE"]] <- subset(d2, d2["hugepages"]==0 & d2["sse"]==1 & d2["popcount_tl"]==0 )
data[["HP"]] <- subset(d2, d2["hugepages"]==1 & d2["sse"]==1 & d2["popcount_tl"]==0 )

form_table <- function(d, order=NA){
	d <- aggregate(d[c('Time','Space')], 
			       by=list(program=d[['program']],test_case=d[['test_case']],
					       hugepages=d[['hugepages']],sse=d[['sse']]),
				   FUN=mean,na.rm=TRUE)
	dByProgram <- split(d, d[["program"]])
	table <- data.frame(dByProgram[[1]]["test_case"])
	names(table) <- c("Text")
	table[["Text"]] <- paste("\\textsc{",table[['Text']],"}") 
	names(table) <- c(" ")
	prog_name <- names(dByProgram)
	if( !is.na(order) ){
		prog_name <- order
	}
	for( prog in prog_name ){
		sel <- dByProgram[[prog]]
		table <- cbind(table, " "=rep("", length(sel["Time"])))
		table <- cbind(table, round(sel["Time"],3))
		table <- cbind(table, sel["Space"]*100)
	}
	unitrow <- paste0(c("", rep(c("&","&($\\mu s$)","&(\\%)"),length(prog_name)), "\\\\[1ex]"),collapse="")
	list("table" = table, "names" = prog_name, "unitrow"=unitrow)
}


make_latex_header <- function(names){
	x <- paste("&&\\multicolumn{2}{c}{", names,"}")
	x <- paste(x, collapse=" ")
	clines=""
	for(i in 1:length(names)){
		clines <- paste(clines,"\\cmidrule{",3*i,"-",3*i+1,"}",sep="")
	}
	y <- paste("\\toprule",x, "\\\\",clines,"\n")	
	gsub("_","\\\\_",y)
}

print_latex <- function( table, names, unitrow ){
	ali <- c("l","r", rep(c("@{\\hspace{1ex}}l","c","c"), (ncol(table)-1)/3) ) 
	dig <- c(0, 0, rep(c(0,3,0),(ncol(table)-1)/3 ))	
	print( xtable( table, align=ali, digits=dig ), 
		   type="latex", hline.after=c(),  # TODO replace by bottomrule
		   floating = F, # don't use table environment
		   add.to.row=list(pos=list(-1,0,nrow(table)), 
				           command=c(make_latex_header(names),unitrow,"\\bottomrule")),
		   sanitize.rownames.function = identity, 
		   sanitize.text.function = identity,
		   include.rownames = FALSE
		 )
}

generate_table <- function(file, data){
	sink(file)
	if ( nrow(data) > 0 ){
		x <- form_table(data, index_ids)
		print_latex(x[["table"]], index_id2latex_name(index_ids), x[["unitrow"]])
	}else{
		cat("\\begin{center}No data for this experiment\\end{center}")
	}
	sink(NULL)
}

sanitize_column <- function(column){
	column <- gsub("_","\\\\_",column)
	column <- gsub("<","{\\\\textless}",column)
	column <- gsub(">","{\\\\textgreater}",column)
	column <- gsub(",",", ",column)
#	column <- gsub(",",",\\\\hspace{0px}",column)
}

for ( feature in names(data) ){
	generate_table(paste0("tbl-count-",feature,".tex"), data[[feature]])
}

sink("tbl-index-info.tex")
info_table <- index_info[c(3,2)]
info_table[[2]] <- sanitize_column(info_table[[2]])

info_table <- cbind(info_table[1], " " = rep("",nrow(info_table)), info_table[2])
ali <- c("l","l","l","p{8cm}")
cat("\\renewcommand{\\arraystretch}{1.3}")
print(xtable(info_table, align=ali), type="latex", hline.after=c(), floating = F, 
	  include.colnames=FALSE, 
	  include.rownames=FALSE, 
	  sanitize.text.function=identity,
	  add.to.row=list(pos=list(-1,nrow(info_table)), 
          			  command=c("\\toprule Identifier&&sdsl class\\\\\\cmidrule{1-1}\\cmidrule{3-3}","\\bottomrule"))
	 ) 
sink(NULL)

