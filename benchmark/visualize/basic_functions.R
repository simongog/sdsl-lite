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
