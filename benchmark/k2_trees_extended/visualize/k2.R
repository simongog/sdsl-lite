require(tikzDevice)
library(gdata)
library(xtable)

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
	mtext("Size relative to original file size (binary adjacency list)", side=1, line=2)

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
plot_time_figure <-function(data,heading,ylab=T,time="ms",xmax=max(data)){
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

	mtext(paste("Time ", time), side=1, line=2)

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

    xmax<-max(data[c('adj_time', 'neighbors_time','reverse_neighbors_time', 'adj_time_comp', 'neighbors_time_comp','reverse_neighbors_time_comp')])
#a <-interleave(data['adj_time'],data['adj_time_comp'])
    #neighbors <-interleave(data['neighbors_time'], data['neighbors_time_comp'])
	#reverse_neighbors <- interleave(data['reverse_neighbors_time'], data['reverse_neighbors_time_comp']);
    #rownames(a)<-interleave(id, paste(id, "_comp"))
    #rownames(neighbors)<-interleave(id, paste(id, "_comp"))
	#rownames(reverse_neighbors)<-interleave(id, paste(id, "_comp"))

a <- data['adj_time']
neighbors <-data['neighbors_time']
reverse_neighbors <-data['reverse_neighbors_time']
rownames(a)<-id
rownames(neighbors)<-id
rownames(reverse_neighbors)<-id


	a_comp <-data['adj_time_comp']
	neighbors_comp <-data['neighbors_time_comp']
	reverse_neighbors_comp <-data['reverse_neighbors_time_comp']
	rownames(a_comp)<-paste(id, "_comp")
	rownames(neighbors_comp)<-paste(id, "_comp")
	rownames(reverse_neighbors_comp)<-paste(id, "_comp")

	time <- "ns"
	if(xmax > 10000){
	    xmax <- xmax/1000
        neighbors <-neighbors/1000
        reverse_neighbors = reverse_neighbors/1000000
	    a <- a/1000
	    neighbors_comp <-neighbors_comp/1000
		reverse_neighbors_comp = reverse_neighbors_comp/1000
		a_comp <- a_comp/1000
		time <- "ms"
	}

	#adj-plot
	#a_complete <- rbind(data.frame(a, index = 1:nrow(a)), data.frame(a_comp, index = 1:nrow(a_comp)))
	#a_complete <- df[order(df$index)]
	plot_time_figure(t(a),"\\tt{adj}", time=time)
    #neighbors-plot
	#neighbors_complete <- rbind(data.frame(neighbors, index = 1:nrow(neighbors)), data.frame(neighbors_comp, index = 1:nrow(neighbors_comp)))
	#neighbors_complete  <- order(neighbors_complete$index)
    plot_time_figure(t(neighbors ),"\\tt{neighbors}", time=time)
    #reverse_neighbors-plot
	#reverse_neighbors_complete <- rbind(data.frame(reverse_neighbors, index = 1:nrow(reverse_neighbors)), data.frame(reverse_neighbors_comp, index = 1:nrow(reverse_neighbors_comp)))
	#reverse_neighbors_complete  <- order(reverse_neighbors_complete$index)
    plot_time_figure(t(reverse_neighbors),"\\tt{reverse_neighbors}", time=time)

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
	plot_time_figure(t(con),"\\tt{construct}",time="sec")

	#construction-size-plot
	tsize<-data[[1,'TC_SIZE']]
	consize <-(data['constructs_space']/tsize)*100
	rownames(consize)<-id

	plot_size_figure(t(consize),"\\tt{construction_space}", ylab=T)

	#size-plot
	tsize<-data[[1,'TC_SIZE']]
	size <-(data['uncompressed_size']/tsize)*100
	rownames(size)<-id
	plot_size_figure(t(size),"\\tt{space}", ylab=T)


	dev.off()
	tex_doc <- paste(tex_doc,"\\begin{figure}[H]
						 \\input{",fig_name,"}
						 \\end{figure}")
	#second page end

	#third page start
	fig_name <- paste("fig-page3-",tc,".tex",sep="")
	open_tikz( fig_name )
	layout(matrix(c(1, 2, 3), 3, 1, byrow=TRUE),
	widths=c(1,1,1), heights=c(1))
	#constructor-plot
	con_comp <-data['compression_time']
	rownames(con)<-id
	plot_time_figure(t(con_comp),"\\tt{compression_time}",time="sec")

	#construction-size-plot
	consize_comp <-(data['compression_space']/tsize)*100
	rownames(consize_comp)<-id

	plot_size_figure(t(consize_comp),"\\tt{compression_space}", ylab=T)

	#size-plot
	size_comp <-(data['compressed_size']/tsize)*100
	rownames(size_comp)<-id
	plot_size_figure(t(size_comp),"\\tt{compressed_size}", ylab=T)

	#third page end
	dev.off()
	tex_doc <- paste(tex_doc,"\\begin{figure}[H]
							 \\input{",fig_name,"}
							 \\end{figure}")

}


#Cool comparison table, do it!
transposed <- maindata
transposed$K2_ID <- NULL
transposed$TC_ID <- NULL
transposed$constructs_time <- paste(transposed$constructs_time," s")
transposed$constructs_space <- paste(formatC(transposed$constructs_space/1024/1024,digits=2,format="f")," Mb")
transposed$constructs_space_vmem <- paste(formatC(transposed$constructs_space_vmem/1024/1024,digits=2,format="f")," Mb")
transposed$construct_sort_duration <- paste(formatC(transposed$construct_sort_duration/1000,digits=2,format="f")," s")
transposed$construct_duration <- paste(formatC(transposed$construct_duration/1000,digits=2,format="f")," s")
transposed$buildvec_duration <- paste(formatC(transposed$buildvec_duration/1000,digits=2,format="f")," s")
transposed$subtree_constructor_duration <- paste(formatC(transposed$subtree_constructor_duration/1000,digits=2,format="f")," s")
transposed$uncompressed_size <- paste(formatC(transposed$uncompressed_size/1024/1024,digits=2,format="f")," Mb")
transposed$adj_time <- paste(transposed$adj_time," ns")
transposed$neighbors_time <- paste(formatC(transposed$neighbors_time/1000,digits=2,format="f")," $\\mu$sec")
transposed$reverse_neighbors_time <- paste(formatC(transposed$reverse_neighbors_time/1000,digits=2,format="f")," $\\mu$sec")

transposed$compression_space <- paste(formatC(transposed$compression_space/1024/1024,digits=2,format="f")," Mb")
transposed$compression_space_vmem <- paste(formatC(transposed$compression_space_vmem/1024/1024,digits=2,format="f")," Mb")
transposed$adj_time_comp <- paste(transposed$adj_time_comp," ns")
transposed$neighbors_time_comp <- paste(formatC(transposed$neighbors_time_comp/1000,digits=2,format="f")," $\\mu$sec")
transposed$reverse_neighbors_time_comp <- paste(formatC(transposed$reverse_neighbors_time_comp/1000,digits=2,format="f")," $\\mu$sec")
transposed$compressed_size <- paste(formatC(transposed$compressed_size/1024/1024,digits=2,format="f")," Mb")
transposed$compression_time <- paste(transposed$compression_time," s")
transposed <- t(transposed)

tex_doc<-paste(tex_doc,print(xtable(transposed), floating = TRUE, floating.environment = "sidewaystable"))

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
