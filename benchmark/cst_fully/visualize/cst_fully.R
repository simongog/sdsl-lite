source("../../basic_functions.R")

idx_config <- readConfig("../index.config", c("IDX_ID","TYPE","LATEX-NAME"))
tc_config <- readConfig("../test_case.config", c("TC_ID","PATH","LATEX-NAME","URL"))

operations = data.frame(id=c("LCA", "LETTER", "SLINK", "CHILD", "DEPTH", "PARENT"),
						caption=c("LCA", "Letter", "SLink", "Child", "SDepth", "Parent"))
rownames(operations)<-operations[["OP_ID"]]
colnames(operations)<-c("OP_ID", "LATEX-NAME")

raw <- data_frame_from_key_value_pairs("../results/all.txt")

sink("fig-space.tex")

cat("\\begin{tabular}{l|*{", nrow(raw), "}{r}}\n", sep="")
cat("\\hline\n")

header = c("")
table1 = c(
    "$\\sigma$",
    "$n/2^{20}$",
    "$|S|$",
    "$\\delta$")
table2 = c(
    "CST (MB)",
    "CSA (MB)")

for(i in 1:nrow(raw)) {
	data<-raw[i,]

    header = paste(header, c(
		tc_config[as.character(data[["TC_ID"]]), "LATEX-NAME"],
		idx_config[as.character(data[["IDX_ID"]]), "LATEX-NAME"]
	), sep=" & ")

	table1 = paste(table1, c(
		data[["ALPHABET_SIZE"]],
		sprintf("%.1f", data[["DATA_SIZE"]] / (1024 * 1024)),
		data[["TREE_SIZE"]],
		data[["SAMPLING_FACTOR"]]
	), sep=" & ")

	table2 = paste(table2, c(
        sprintf("%.2f", data[["CST_SIZE"]] / (1024 * 1024)),
        sprintf("%.2f", data[["CSA_SIZE"]] / (1024 * 1024))
	), sep=" & ")
}

cat(header, "\\hline", sep="\\\\\n")
cat(table1, "\\hline", sep="\\\\\n")
cat(table2, "\\hline", sep="\\\\\n")
cat("\\end{tabular}\n")

sink("fig-time.tex")

cat("\\begin{tabular}{l|*{", nrow(raw), "}{r}}\n", sep="")
cat("\\hline\n")

header = c("Operation", "")
for(i in 1:nrow(raw)) {
	data<-raw[i,]
    header = paste(header, c(
		tc_config[as.character(data[["TC_ID"]]), "LATEX-NAME"],
		idx_config[as.character(data[["IDX_ID"]]), "LATEX-NAME"]
	), sep=" & ")
}
cat(header, "\\hline", sep="\\\\\n")

format_string = "%.2f"

for(op_idx in 1:nrow(operations)) {
	op_name = as.character(operations[op_idx, "OP_ID"])
	table = c(paste(as.character(operations[op_idx, "LATEX-NAME"]), sep="", collapse=""))

	for(i in 1:nrow(raw)) {
		data<-raw[i,]

		table = paste(table, c(
			sprintf(format_string, data[[paste("CST_", op_name, "_U_TIME", sep="", collapse="")]] / 1000)
		), sep=" & ")
	}

	cat(table, "\\\\\n", sep="")
}

cat("\\end{tabular}\n")

sink(NULL)
