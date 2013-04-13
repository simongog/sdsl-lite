require("grDevices") # for color

args <- commandArgs(trailingOnly = TRUE)

source(args[1])
input=args[2]
tikz_output=as.numeric(args[3])
divspace=as.numeric(args[4])
timeunit=args[5]
ylabelladd=args[6]

slashes=""
if ( tikz_output ){ 
	if ( !exists( "tikzDeviceLoaded" ) ){  
		require(tikzDevice) 
		tikzDeviceLoaded = T
	}
	tikz("mm-log.tex", width = 3.2, height = 1.3, standAlone = T)
	par(oma=c(0,0,0,0))
	par(mar=c(2.3,2.4,1.2,0.2))
	par(xpd=F)
	par(mgp=c(1,0.3,0))
	par(xaxs="r")
	par(yaxs="i")
	slashes="\\"
}else{
	pdf("mm-log.pdf", width = 3*3.2, height = 3*1.3)	
}

t <- read.csv(input,sep=";",header=F)
t$V2 <- t$V4
t$V3 <- t$V5
t$V2 <- t$V2/divspace
if ( timeunit == "seconds" ) {
	t$V1 <- t$V1/1000.0
}

maxy=max(t$V2*1)
maxx=max(t$V1)
dx = maxx/500
plot(t$V1, t$V2, type="l",col="gray", ylab=paste("Space ",ylabelladd,sep=""), xlab=paste("Time (",timeunit,")",sep=""), bty="l",
	 ylim=c(0, maxy), cex.axis=0.6, cex.lab=0.7,las=1, 
	 tck=-0.03 # length of ticks
	 )
abline(h=axTicks(2), col="gray", lty="dotted")
polygon(x_for_polygon(t$V1), y_for_polygon(t$V2), col="lightgray", border=NA  )
lines(t$V1, t$V2, type="l",col="gray")

config <- readConfig("mm-log.config",c("ID","start","end","fill-color","stroke-color","NAME"))

par(xpd=NA)
phase_cnt <- 1
for( x in 1:nrow(config) ){
	print(config[x,"start"] )
	t1 <- subset(t, t[3]==config[x,"start"])[[1]]
	t2 <- subset(t, t[3]==config[x,"end"])[[1]]
	interval <- subset(t, t[1] >= t1 & t[1] <= t2)
	xs <- interval[[1]]
	ys <- interval[[2]]
	polygon(x_for_polygon(xs), y_for_polygon(ys), col=config[x,"fill-color"], border=NA  )
	lines(xs, ys,col=config[x,"stroke-color"])
	lines(c(rep(head(xs,1)+dx,2),rep(tail(xs,1)-dx,2)), c(maxy, rep(maxy*1.02 ,2), maxy) )
	phase_str = paste(slashes,"#",phase_cnt,sep="")
	if ( tail(xs,1)-head(xs,1) > 2.5*strwidth(phase_cnt, cex=0.6 ) ){
		text(head(xs,1), maxy*1.07, phase_str ,adj=0, cex=0.6)
	}
	phase_cnt <- phase_cnt + 1
}



dev.off()
