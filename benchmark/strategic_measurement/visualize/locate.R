source("../../basic_functions.R")

tc_config <- readConfig("../test_case.config",c("TC_ID","PATH","LATEX-NAME","URL"))

sink("locate.tex")
cat("\\documentclass{article}\
\\usepackage{booktabs}\
\
\\begin{document}")


raws <- data_frame_from_key_value_pairs( "../results/locate.txt")

for ( tc_id in tc_config[["TC_ID"]] ){
tc_name <- tc_config[tc_id,"LATEX-NAME"]

raw <- subset(raws, raws[["TC_ID"]]==tc_id)
raw[["k"]] <- (raw[["PAT_MIN_OCC"]]+raw[["PAT_MAX_OCC"]])/2
raw[["m"]] <- raw[["PAT_M"]]
raw[["time"]] <- raw[["Locate_time_in_secs"]]*1000
raw2 <- raw
raw <- aggregate(raw2,by=list(raw2[["m"]],raw2[["k"]]),FUN=mean,na.rm=T)
raw[["sd_time"]] <- aggregate(raw2,by=list(raw2[["m"]],raw2[["k"]]),FUN=sd,na.rm=T)[["time"]]


raw[["avg_time"]] <- sprintf("$%3.2f$",round(raw[["time"]],2))
raw[["avg_time"]] <- gsub(" ","\\\\hphantom{0}",raw[["avg_time"]])
raw[raw[["file_pat_cnt"]]<10000,"avg_time"] <- "\\multicolumn{1}{c}{-}"
raw[["sd_time2"]] <- sprintf("$%8.2f$",round(raw[["sd_time"]],2))
raw[["sd_time2"]] <- gsub(" ","\\\\hphantom{0}",raw[["sd_time2"]])
raw[raw[["file_pat_cnt"]]<10000,"sd_time2"] <- "\\multicolumn{1}{c}{-}"
raw <- cbind( raw, D = rep("",nrow(raw))) # add space column
raw <- raw[order(raw[["m"]],raw[["k"]]),]
data <- split(raw[c("D","avg_time","sd_time2")], raw[["m"]])
tab <- do.call("cbind", data)

ks <- subset(raw,raw[["m"]]==8)[["k"]]
kcol <- gsub("$k=0$","uncontrolled",paste("$k=",ks,"$",sep=""),fixed=T)

tab <- cbind(kcol, tab)

mypaste<-function(...){paste(...,sep="&")}
cat("\\begin{table}\n\\centering\n")
cat("\\begin{tabular}{l*{",length(data),"}{@{\\quad}l@{}r@{$\\,\\pm\\,$}c}}\n")
cat("\\toprule\n")
cat("Answers $k$",paste("&&\\multicolumn{2}{c}{$m=",names(data),"$}",sep=""),"\\\\")
cat("\\cmidrule{1-",length(data)*3+1,"}\n",sep="")
cat(paste(do.call("mypaste", tab),collapse="\\\\\n"))
cat("\\\\\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\caption{Time in milliseconds for locate on different pattern sets for test case ",tc_name,".}\n")
cat("\\end{table}\n\n")

}

cat("\\end{document}")
sink(NULL)
