require(tikzDevice)
source("../../basic_functions.R")

tex_file = "k2.tex"

tc_config <- readConfig("../test_case.config",c("TC_ID","PATH","LATEX_NAME","URL"))

open_tikz <- function( file_name ){
    tikz(file_name, width = 5.5, height = 7.5 , standAlone = F , sanitize = TRUE)
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
	mtext("Size relative to original file size (arc file)", side=1, line=2)

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
tex_doc <- paste(readLines("k2-header.tex"),collapse="\n")

tex_doc<-paste(tex_doc,"\\section{Result of the K2 Tree benchmark}")


maindata <- data_frame_from_key_value_pairs( "../results/all.txt" )

#create two pages for each test case
#for(tc in tc_config[['TC_ID']]){
for(tc in unique(maindata$TC_ID)){

	data<-maindata[maindata$TC_ID==tc,]
	id <-data[['K2_TEX_NAME']]

	#first page start
	fig_name <- paste("fig-page1-",tc,".tex",sep="")
	tex_doc<-paste(tex_doc,"\\subsection{Test case: {\\sc ",data[['TC_TEX_NAME']],"}}")

	open_tikz( fig_name )

	layout(matrix(c(1,2,3), nrow=3, ncol=1, byrow = TRUE),
	   widths=c(1,1,1), heights=c(1))

    xmax<-max(data[c('adj_time', 'neighbors_time','reverse_neighbors_time')])
	a <-data['adj_time']
    neighbors <-data['neighbors_time']
    reverse_neighbors <-data['reverse_neighbors_time']
    rownames(a)<-id
    rownames(neighbors)<-id
    rownames(reverse_neighbors)<-id
	if(xmax > 10000){
	    xmax <- xmax/1000000
        neighbors <-neighbors/1000000
        reverse_neighbors = reverse_neighbors/1000000
	    a <- a/1000000
	    #adj-plot
	    plot_time_figure(t(a),"\\tt{adj}", xlab=F)
        #neighbors-plot
        plot_time_figure(t(neighbors),"\\tt{neighbors}", xlab=F, xmax=xmax)
        #reverse_neighbors-plot
        plot_time_figure(t(reverse_neighbors),"\\tt{reverse_neighbors}",constructor=T, xlab=F, xmax=xmax)
    }
    else {
	    #adj-plot
	    plot_time_figure(t(a),"\\tt{adj}", xlab=F)
        #neighbors-plot
        plot_time_figure(t(neighbors),"\\tt{neighbors}", xlab=F, xmax=xmax)
        #reverse_neighbors-plot
        plot_time_figure(t(reverse_neighbors),"\\tt{reverse_neighbors}", xmax=xmax)

    }

    old<-par()
    dev.off()
    tex_doc <- paste(tex_doc,"\\begin{figure}[H]
                    \\input{",fig_name,"}
                    \\end{figure}")
	#first page end

	#second page start
	fig_name <- paste("fig-page2-",tc,".tex",sep="")
	open_tikz( fig_name )

	layout(matrix(c(1, 2, 3), 3, 1, byrow=TRUE),
	   widths=c(1,1,1), heights=c(1))

	#constructor-plot
	con <-data['constructs_time']
	rownames(con)<-id
	plot_time_figure(t(con),"\\tt{construct}",xlab=F,constructor=T)

	#construction-size-plot
	tsize<-data[[1,'TC_SIZE']]
	consize <-(data['constructs_space']/tsize)*100
	rownames(consize)<-id

	plot_size_figure(t(consize),"\\tt{construction space}", ylab=T)

	#size-plot
	tsize<-data[[1,'TC_SIZE']]
	size <-(data['k2_size']/tsize)*100
	rownames(size)<-id
	plot_size_figure(t(size),"\\tt{space}", ylab=T)

	dev.off()
	tex_doc <- paste(tex_doc,"\\begin{figure}[H]
					 \\input{",fig_name,"}
					 \\end{figure}")
	#second page end
}

#type identification table
tex_doc<-paste(tex_doc,"\\begin{table}[b]
						\\centering",
						typeInfoTable("../k2tree.config",data[['K2_ID']], 1, 3, 2),
						"\\caption{K2 tree identifier and corresponding sdsl-type.}
						\\end{table}")

#read footer+end
tex_doc <- paste(tex_doc, readLines("k2-footer.tex"),collapse="\n")
sink(tex_file)
cat(tex_doc)
sink(NULL)
