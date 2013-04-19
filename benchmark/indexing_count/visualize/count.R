library(xtable) # if not installed call install.packages("xtable")
library(plyr)

source("../../basic_functions.R")

# Load index information
idx_config <- readConfig("../index.config",c("IDX_ID","SDSL_TYPE","LATEX-NAME"))
tc_config <- readConfig("../test_case.config",c("TC_ID","PATH","LATEX-NAME","URL"))
compile_config <- readConfig("../compile_options.config",c("COMPILE_ID","OPTIONS"))

# Create data frame which maps test cases names to their index in the list
tc_ord <- data.frame("ord"=seq(1,nrow(tc_config)),"LATEX-NAME")
rownames(tc_ord) <- tc_config[["TC_ID"]]

# Load report information

config <- readConfig("index-filter.config",c("IDX_ID"))

# Load data
raw <- data_frame_from_key_value_pairs( "../results/all.txt" )
#
# Filer indexes
raw               <- raw[raw[["IDX_ID"]]%in%config[["IDX_ID"]],] 
raw[["IDX_ID"]] <- factor(raw[["IDX_ID"]])
# Normalize data
raw[["Time"]]     <- raw[["Count_time_in_milli_sec"]]*1000
raw[["Time"]]     <- raw[["Time"]]/(raw[["pat_cnt"]]*raw[["pat_length"]])
raw[["Time"]]     <- round(raw[["Time"]], 3)
raw[["Space"]]    <- raw[["Index_size_in_bytes"]]/raw[["text_size"]]
raw[["Space"]]    <- round(raw[["Space"]],2)

raw <- raw[c("TC_ID", "Space", "Time","COMPILE_ID","IDX_ID")]
raw <- raw[order(raw[["TC_ID"]]),]

data <- split(raw, raw[["COMPILE_ID"]])



form_table <- function(d, order=NA){
# calculate the mean time per IDX_ID,TC_ID
    d <- aggregate(d[c('Time','Space')], 
                   by=list(IDX_ID=d[['IDX_ID']],
                           TC_ID=d[['TC_ID']]),
                   FUN=mean,na.rm=TRUE)
	d <- d[ order(tc_ord[as.character(d[["TC_ID"]]),"ord"]), ]
    dd <- split(d, d[["IDX_ID"]])
    table <- data.frame(dd[[1]]["TC_ID"], stringsAsFactors=F)
    names(table) <- c("TC_ID")
    table[["TC_ID"]] <- paste("\\textsc{", tc_config[as.character(table[['TC_ID']]), "LATEX-NAME"],"}") 
    names(table) <- c(" ")
    prog_name <- names(dd)
    if( !is.na(order) ){
        prog_name <- order
    }
    for( prog in prog_name ){
        sel <- dd[[prog]]
        table <- cbind(table, " "=rep("", length(sel["Time"])))
        table <- cbind(table, round(sel["Time"],3))
        table <- cbind(table, sel["Space"]*100)
    }
    unitrow <- paste(c("", rep(c("&","&($\\mu s$)","&(\\%)"),length(prog_name)), "\\\\[1ex]"),collapse="",sep='')
    list("table" = table, "names" = table[["TC_ID"]], "unitrow"=unitrow)
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
        x <- form_table(data, config[["IDX_ID"]])
        print_latex(x[["table"]], idx_config[as.character(config[["IDX_ID"]]), "LATEX-NAME"], x[["unitrow"]])
    }else{
        cat("\\begin{center}No data for this experiment\\end{center}")
    }
    sink(NULL)
}

for ( compile_id in names(data) ){
    generate_table(paste("tbl-count-",compile_id,".tex",sep=''), data[[compile_id]])
}

sink("count.tex")
cat(paste(readLines("count-header.tex"),"\n",sep=""))
for ( compile_id in names(data) ){
    cat("
    \\begin{table}
    \\centering
        \\input{tbl-count-",compile_id,".tex}
    \\caption{Time in $\\mu$sec per pattern symbol in a count query.
             Index space as fraction of original file size.
             Compile options: 
             \\texttt{",gsub("_","\\\\_",compile_config[compile_id, "OPTIONS"]),"}.
             \\label{tbl-count-",compile_id,"}}
    \\end{table}",sep="")
}

# output table containing type information
cat("\\begin{table}
    \\centering",
typeInfoTable("../index.config", config[["IDX_ID"]], 1, 3, 2)
    ,"\\caption{Index identifier and corresponding sdsl-type.}
    \\end{table}
    ")

cat(paste(readLines("count-footer.tex"),"\n",sep=""))
sink(NULL)
