# Read a file called file_name and create a data frame the following way
# (1) Parse all the lines of the form
# '# key = value'
# (2) Each unique key gets a column
data_frame_from_key_value_pairs <- function(file_name){
    lines <- readLines(file_name)
    lines <- lines[grep("^#.*=.*",lines)]
    d <- gsub("^#","",gsub("[[:space:]]","",unlist(strsplit(lines,split="="))))
    keys <- unique(d[seq(1,length(d),2)])
    keynr <- length(keys)
    dd <- d[seq(2,length(d),2)]
    dim(dd) <- c( keynr, length(dd)/keynr )
    data <- data.frame(t(dd))    
    names(data) <- keys    

    for (col in keys){
        t <- as.character(data[[col]])
        suppressWarnings( tt <- as.numeric(t) )
        if ( length( tt[is.na(tt)] ) == 0 ){ # if there are not NA in tt
            data[[col]] <- tt 
        }
    }     
    data
}

# Takes a vector v=(v1,v2,v3,....)
# and returns a vector which repeats
# each element x times. So for two we get
# (v1,v1,v2,v2,v3,v3....)    
expand_vec <- function( v, x ){
    v <- rep(v,x)
    dim(v) <- c(length(v)/x,x)
    v <- t(v)
    dim(v) <- c(1,nrow(v)*ncol(v))
    v
}

# Takes a vector v=(v1,v2,v3,....)
# and returns a vector which appends x-1
# NA after each value. So for x=2 we get
# (v1,NA,v2,NA,v3,NA....)    
expand_vec_by_NA <- function( v, x ){
    v <- c(v, rep(NA,length(v)*(x-1)))
    dim(v) <- c(length(v)/x,x)
    v <- t(v)
    dim(v) <- c(1,nrow(v)*ncol(v))
    v
}

format_str_fixed_width <- function(x, width=4){
    sx <- as.character(x)
    if ( nchar(sx) < width ){
        for (i in 1:(width-nchar(sx))){
            sx <- paste("\\D",sx, sep="")
        }
    }
    sx
}

# Check if package is installed
# found at: http://r.789695.n4.nabble.com/test-if-a-package-is-installed-td1750671.html#a1750674
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 


sanitize_column <- function(column){
    column <- gsub("_","\\\\_",column)
    column <- gsub(" >",">",column)
    column <- gsub("<","{\\\\textless}",column)
    column <- gsub(">","{\\\\textgreater}",column)
    column <- gsub(",",", ",column)
}

# tranforms a vector of index ids to a vector which contains the
# corresponding latex names.
# Note: each id should only appear once in the input vector
mapids <- function(ids, mapit){
    as.character(unlist(mapit[ids]))
}

id2latex <- function(config_file, latexcol, idcol=1){
  index_info <- read.csv(config_file, sep=";",header=F, comment.char="#")
  res <- data.frame( t(as.character(index_info[[latexcol]])), stringsAsFactors=F )
  names(res) <- as.character(index_info[[idcol]]) 
  res 
}

idAndValue <- function(config_file, valuecol, idcol=1){
  res <- read.csv(config_file, sep=";",header=F, comment.char="#")
  res[c(idcol, valuecol)]
}

readConfig <- function(config_file, mycolnames){
  config <- read.csv(config_file, sep=";",header=F, comment.char="#",stringsAsFactors=F)
  rownames(config) <- config[[1]]
  colnames(config) <- mycolnames
  config
}

# Creates a LaTeX table containing index names and sdsl type
# config_file   The index.config storing the type information
# index_ids     Filter the index.config entires with this index ids
# id_col        Column `id_col` contains the IDs
# name_col      Column `name_col` contains the latex names
# type_col      Column `type_col` contains the type
typeInfoTable <- function(config_file, index_ids, id_col=1, name_col=3, type_col=2){
    x <- read.csv(config_file, sep=";", header=F, comment.char="#",stringsAsFactors=F)
    rownames(x) <- x[[id_col]]
    x <- x[index_ids,] # filter
    sdsl_type <- sanitize_column(x[[type_col]])
    sdsl_name <- x[[name_col]]
    res <- "
        \\renewcommand{\\arraystretch}{1.3}
        \\begin{tabular}{@{}llp{10cm}@{}}
          \\toprule
          Identifier&&sdsl type\\\\ \\cmidrule{1-1}\\cmidrule{3-3}"
    res <- paste(res, paste(sdsl_name,"&&\\footnotesize\\RaggedRight\\texttt{",sdsl_type,"}\\\\",sep="",collapse=" "))
    res <- paste(res,"
        \\bottomrule
        \\end{tabular}")
}

# returns x concatenated with x reversed 
x_for_polygon <- function(x){
  c( x, rev(x) )
}

# return y concatenated with rep(0, length(y))
y_for_polygon <- function(y){
  c( y, rep(0, length(y)) )
}

# ncols Number of columns in the figure
# nrows Number of rows in the figure
multi_figure_style <- function(nrows, ncols){
  par(mfrow=c(nrows, ncols))

  par(las=1) # axis labels always horizontal
  par(yaxs="i") # don't add +- 4% to the yaxis
  par(xaxs="i") # don't add +- 4% to the xaxis


  # distance (x1,x2,x3) of axis parts from the axis. x1=axis labels or titles
  # x2=tick marks, x3=tick marks symbol
  par(mgp=c(2,0.5,0)) 

  # length of tick mark as a fraction of the height of a line of text, default=-0.5
  par(tcl=-0.2) 

  par(oma=c(2.5,2.7,0,0.2)) # outer margin (bottom,left,top,right)
  par(mar=c(1,1,1.5,0.5)) # inner margin (bottom,left,top,right)

}

# Draw the heading of diagrams
# text  Text which should be displayed in the heading
draw_figure_heading <- function(text){
    # scale Y
    SY <- function(val){ if( par("ylog") ){ 10^val } else { val } }
    SX <- function(val){ if( par("xlog") ){ 10^val } else { val } }

    rect(xleft=SX(par("usr")[1]), xright=SX(par("usr")[2]), 
         ybottom=SY(par("usr")[4]), ytop=SY(par("usr")[4]*1.1) ,xpd=NA,
         col="grey80", border="grey80" )
    text(labels=text,y=SY(par("usr")[4]*1.02), adj=c(0.5, 0),x=SX((par("usr")[1]+par("usr")[2])/2),xpd=NA,cex=1.4)
}

print_info <- function(){
	Sys.info("release")
	Sys.info("sysname")
	Sys.info("version")
	Sys.info("nodename")
	Sys.info("machine")
	Sys.info("login")
	Sys.info("user")
	Sys.info("effective_user")
}
