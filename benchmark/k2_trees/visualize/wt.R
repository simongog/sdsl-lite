require(tikzDevice)
source("../../basic_functions.R")

tex_file = "wt.tex"

tc_config <- readConfig("../test_case.config",c("TC_ID","PATH","LATEX_NAME","URL","TC_TYPE"))

open_tikz <- function( file_name ){
    tikz(file_name, width = 5.5, height = 7.5 , standAlone = F) 
}

x_for_bar<-function(value){
	c(0,0,value,value)
}

y_for_bar<-function(offset){
	c(offset,offset+0.4,offset+0.4,offset)	
}

#Method which plots the size figure
plot_size_figure <-function(data,heading,ylab=F){

	#set margin
	par(mar=c(3,2,3,0))
	if(ylab){
		par(mar=c(3,10,3,0))
	}

	plot(c(),c(),ylim=c(0,(length(data)*0.5)+0.2),xlim=c(0,max(101,(max(data)+1))),xlab="",ylab="",xaxt="n",yaxt="n")

	#label y-axis
	if(ylab){
		axis( 2, at =seq(0.3,(length(data)*0.5)+0.2,0.5), label=colnames(data),las=1)
	}
	#label x-axis
	axis(1)
	mtext("Size relative to original file size", side=1, line=2)

	#draw bars
	offset=0.1
	for(time in data){
		polygon( x_for_bar(time),y_for_bar(offset), border=NA, col="grey")
		offset=offset+0.5
	}

    #abline(v=c(axis(1)/2,max(axis(1)/2)+axis(1)/2), col="gray")
    abline(v=c(axis(1),axis(1)+(axis(1)[2]-axis(1)[1])/2),col="gray")
	abline(v=100, col="red")
	draw_figure_heading(heading)
}

#Method which plots the a time figure
plot_time_figure <-function(data,heading,ylab=T,xlab=T,constructor=F,xmax=max(data)){
	#set margin
	par(mar=c(3,2,2,0))
	if(ylab){
		par(mar=c(3,10,2,0))
	}

	plot(c(),c(),ylim=c(0,(length(data)*0.5)+0.2),xlim=c(0,(xmax*1.02)),xlab="",ylab="",xaxt="n",yaxt="n")

	#label y-axis
	if(ylab){
		axis( 2, at =seq(0.3,(length(data)*0.5)+0.2,0.5), label=colnames(data),las=1)		
	}	
	#label x-axis
	axis(1)
    abline(v=c(axis(1),axis(1)+(axis(1)[2]-axis(1)[1])/2),col="gray")
	if(xlab){
		mtext("Time in microseconds", side=1, line=2)
	}
	if(constructor){
		mtext("Time in seconds", side=1, line=2)
	}

	#draw bars
	offset=0.1
	for(time in data){
		polygon( x_for_bar(time),y_for_bar(offset), border=NA, col="grey")
		offset=offset+0.5
	}

	draw_figure_heading(heading)
}


#read header
tex_doc <- paste(readLines("wt-header.tex"),collapse="\n")

tex_doc<-paste(tex_doc,"\\section{Result of the Wavelet Tree benchmark}")


maindata <- data_frame_from_key_value_pairs( "../results/all.txt" )

#create two pages for each test case
#for(tc in tc_config[['TC_ID']]){
for(tc in unique(maindata$TC_ID)){

	data<-maindata[maindata$TC_ID==tc,]
	id <-data[['WT_TEX_NAME']]
	xmax<-max(data[c('access_time','rank_time','select_time','inverse_select_time','lex_count_time','lex_smaller_count_time')])

	#first page start 
	fig_name <- paste("fig-page1-",tc,".tex",sep="")
	tex_doc<-paste(tex_doc,"\\subsection{Test case: {\\sc ",data[['TC_TEX_NAME']],"}}")

	open_tikz( fig_name )

	layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE),
	   widths=c(1.35,1), heights=c(1,1,1))

	#access-plot
	a <-data['access_time']
	rownames(a)<-id
	plot_time_figure(t(a),"\\tt{access}",xlab=F,xmax=xmax)

	#rank-plot
	rank <-data['rank_time']
	rownames(rank)<-id
	plot_time_figure(t(rank),"\\tt{rank}",ylab=F,xlab=F,xmax=xmax)

	#select-plot
	s <-data['select_time']
	rownames(s)<-id
	plot_time_figure(t(s),"\\tt{select}",xlab=F,xmax=xmax)

	#inverse-select-plot
	is <-data['inverse_select_time']
	rownames(is)<-id
	plot_time_figure(t(is),"\\tt{inverse\\_select}",xlab=F,ylab=F,xmax=xmax)

	#lex-count-plot
	lc <-data['lex_count_time']
	rownames(lc)<-id
	plot_time_figure(t(lc),"\\tt{lex\\_count}",xmax=xmax)

	#lex-smaller-count-plot
	lsc <-data['lex_smaller_count_time']
	rownames(lsc)<-id
	plot_time_figure(t(lsc),"\\tt{lex\\_smaller\\_count}",ylab=F,xmax=xmax)
	
	old<-par()
	dev.off()
	tex_doc <- paste(tex_doc,"\\begin{figure}[H]
					 \\input{",fig_name,"}
					 \\end{figure}")
	#first page end

	#second page start
	fig_name <- paste("fig-page2-",tc,".tex",sep="")
	open_tikz( fig_name )

	layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE),
	   widths=c(1.35,1), heights=c(1,1,1))

	#interval-symbols-plot
	ivs <-data['interval_symbols_time']
	rownames(ivs)<-id
	plot_time_figure(t(ivs),"\\tt{interval\\_symbols}",xmax=max(xmax,max(ivs)))

	#constructor-plot
	con <-data['constructs_time']
	rownames(con)<-id
	plot_time_figure(t(con),"\\tt{construct}",ylab=F,xlab=F,constructor=T)

	#construction-size-plot
	tsize<-data[[1,'TC_SIZE']]
	consize <-(data['constructs_space']/tsize)*100	
	rownames(consize)<-id

	plot_size_figure(t(consize),"\\tt{construction space}",ylab=T)

	#size-plot
	tsize<-data[[1,'TC_SIZE']]
	size <-(data['wt_size']/tsize)*100
	rownames(size)<-id

	plot_size_figure(t(size),"\\tt{space}")

	dev.off()
	tex_doc <- paste(tex_doc,"\\begin{figure}[H]	
					 \\input{",fig_name,"}
					 \\end{figure}")
	#second page end
}

#type identification table
tex_doc<-paste(tex_doc,"\\begin{table}[b]
						\\centering",
						typeInfoTable("../wt.config",data[['WT_ID']], 1, 3, 2),
						"\\caption{Wavelet tree identifier and corresponding sdsl-type.}
						\\end{table}")

#read footer+end
tex_doc <- paste(tex_doc, readLines("wt-footer.tex"),collapse="\n")
sink(tex_file)
cat(tex_doc)
sink(NULL)
