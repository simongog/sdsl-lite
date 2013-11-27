library(xtable)
source("../../basic_functions.R")

tex_file = "lcp.tex"

tc_config <- readConfig("../test_case.config",c("TC_ID","PATH","LATEX_NAME","URL"))
lcp_config <- readConfig("../lcp.config",c("LCP_ID","LCP_TYPE","LATEX_NAME","BWT"))


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

#read header
sink(tex_file)
cat(paste(readLines("lcp-header.tex"),collapse="\n"))

maindata <- data_frame_from_key_value_pairs( "../results/all.txt" )

names<-c("SA","BWT","LCP","OVERALL")
unitrow <- paste(c("", rep(c("&&Time", "&Space"), length(names)), "\\\\","", rep(c("&&(sec)", "&(\\%)"), length(names)), "\\\\[1ex]"), collapse="", sep='')

# create a table for each test case
for(i in 1:nrow(maindata)){

	data<-maindata[i,]
	row<-nrow(lcp_config)
	size<-data[["TC_SIZE"]]
	table<-data.frame(EMPTY=character(row),SATIME=character(row),SASPACE=character(row),EMPTY2=character(row),BWTTIME=character(row),BWTSPACE=character(row),EMPTY3=character(row),LCPTIME=character(row),LCPSPACE=character(row),EMPTY4=character(row),OVERALLTIME=character(row),OVERALLSPACE=character(row),stringsAsFactors=FALSE)

	# gather data
	for(l in 1:row){
		table[l,]["SATIME"]<-sprintf("%.2f",data[["SA_TIME"]])
		table[l,]["SASPACE"]<-round(data[["SA_MMPEAK"]]*100/size, digits=0)


		if(lcp_config[["BWT"]][l]){
			table[l,]["BWTTIME"]<-sprintf("%.2f",data[["BWT_TIME"]])
			table[l,]["BWTSPACE"]<-round(data[["BWT_MMPEAK"]]*100/size, digits=0)
			table[l,]["OVERALLTIME"]<-sprintf("%.2f",data[["SA_TIME"]]+data[["BWT_TIME"]]+data[[paste(lcp_config[["LCP_ID"]][l],"_TIME",sep="")]])
			table[l,]["OVERALLSPACE"]<-round(max(data[["SA_MMPEAK"]],data[["BWT_MMPEAK"]],data[[paste(lcp_config[["LCP_ID"]][l],"_MMPEAK",sep="")]])*100/size, digits=0)
		}
		else{
			table[l,]["BWTTIME"]<-"-"
			table[l,]["BWTSPACE"]<-"-"
			table[l,]["OVERALLTIME"]<-sprintf("%.2f",data[["SA_TIME"]]+data[[paste(lcp_config[["LCP_ID"]][l],"_TIME",sep="")]])
			table[l,]["OVERALLSPACE"]<-round(max(data[["SA_MMPEAK"]],data[[paste(lcp_config[["LCP_ID"]][l],"_MMPEAK",sep="")]])*100/size, digits=0)
		}

		table[l,]["LCPTIME"]<-sprintf("%.2f",data[[paste(lcp_config[["LCP_ID"]][l],"_TIME",sep="")]])
		table[l,]["LCPSPACE"]<-round(data[[paste(lcp_config[["LCP_ID"]][l],"_MMPEAK",sep="")]]*100/size, digits=0)
	}

	row.names(table)<-lcp_config[["LATEX_NAME"]]

	# convert and print table
	ali <- c("l", rep(c("@{\\hspace{1ex}}l","c","c"), (ncol(table))/3) )
	dig <- c(0,  rep(c(0,3,0),(ncol(table))/3 ))

	print(	xtable(table, align=ali, digits=dig,
			caption = paste("Results for ",as.character(data[["TC_TEX_NAME"]])," (size: ",round(size/(1024^2), digits=3),"MB). Runtime in seconds. Space is the peak memory usage (including input and output) as fraction of original file size.")),
           	add.to.row=list(pos=list(-1,0,nrow(table)), command=c(make_latex_header(names),unitrow,"\\bottomrule")),
		hline.after=c(),
           	sanitize.rownames.function = identity,
		include.colnames = FALSE
	)
}

cat(paste(readLines("lcp-footer.tex"),collapse="\n"))
sink(NULL)
