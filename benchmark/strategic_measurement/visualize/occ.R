source("../../basic_functions.R")

tc_config <- readConfig("../test_case.config",c("TC_ID","PATH","LATEX-NAME","URL"))

raw <- data_frame_from_key_value_pairs( "../results/occ.txt")
raw[["avg_occ"]] <- round(raw[["avg_occ"]],2)
raw <- cbind( raw, D = rep("",nrow(raw))) # add space column

data <- split(raw[c("D","avg_occ","med_occ")], raw[["m"]])

tab <- do.call("cbind", data)
tab <- cbind( name=tc_config[as.character(raw[rownames(tab),"TC_ID"]),"LATEX-NAME"]  , tab)

sink("tbl-occ.tex")
mypaste<-function(...){paste(...,sep="&")}
cat("\\begin{tabular}{r*{",length(data),"}{lr@{\\ }r}}\n")
cat("\\toprule\n")
cat(paste("&&\\multicolumn{2}{c}{$m=",names(data),"$}",sep=""),"\\\\")
for(i in 1:length(data)){cat("\\cmidrule{",3*i,"-",3*i+1,"}",sep="")}
cat(paste("&&",rep("\\multicolumn{1}{c}{avg} & \\multicolumn{1}{c}{med}",length(data)),sep=""),"\\\\[1ex]")
cat(paste(do.call("mypaste", tab),collapse="\\\\\n"))
cat("\\\\\\bottomrule\n")
cat("\\end{tabular}\n")
sink(NULL)
