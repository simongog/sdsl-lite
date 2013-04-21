source("../../basic_functions.R")

tc_config <- readConfig("../test_case.config",c("TC_ID","PATH","LATEX-NAME","URL"))

sink("cand.tex")
cat("\\documentclass{article}\
\\usepackage{booktabs}\
\
\\begin{document}")


raws <- data_frame_from_key_value_pairs( "../results/cand.txt")

for ( tc_id in tc_config[["TC_ID"]] ){
tc_name <- tc_config[tc_id,"LATEX-NAME"]

raw <- subset(raws, raws[["TC_ID"]]==tc_id)
raw[["k"]] <- (raw[["PAT_MIN_OCC"]]+raw[["PAT_MAX_OCC"]])/2
raw <- subset(raw, raw["k"]>0)
raw[["m"]] <- raw[["PAT_M"]]
raw[["cand"]] <- sprintf("$%10d$",raw[["CANDIDATES"]])
raw[["cand"]] <- gsub(" ","\\\\hphantom{0}",raw[["cand"]])
#raw[raw[["CANDIDATES"]]<0,"cand"] <- "\\multicolumn{1}{c}{-}"
raw <- cbind( raw, D = rep("",nrow(raw))) # add space column
raw <- raw[order(raw[["m"]],raw[["k"]]),]
data <- split(raw[c("D","cand")], raw[["m"]])
tab <- do.call("cbind", data)

ks <- subset(raw,raw[["m"]]==8)[["k"]]
kcol <- gsub("$k=0$","uncontrolled",paste("$k=",ks,"$",sep=""),fixed=T)

tab <- cbind(kcol, tab)

mypaste<-function(...){paste(...,sep="&")}
cat("\\begin{table}\n\\centering\n")
cat("\\begin{tabular}{l*{",length(data),"}{lc}}\n")
cat("\\toprule\n")
cat("Answers $k$",paste("&&\\multicolumn{1}{c}{$m=",names(data),"$}",sep=""),"\\\\")
cat("\\cmidrule{1-",length(data)*2+1,"}\n",sep="")
cat(paste(do.call("mypaste", tab),collapse="\\\\\n"))
cat("\\\\\\bottomrule\n")
cat("\\end{tabular}\n")
cat("\\caption{Number of pattern candidates in ",tc_name,".}\n")
cat("\\end{table}\n\n")

}

cat("\\end{document}")
sink(NULL)
