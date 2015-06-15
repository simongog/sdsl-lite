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
    "$|T|/2^{20}$")
table2 = c()
for(idx_idx in 1:nrow(idx_config)) {
	table2 = c(table2, paste(
		idx_config[as.character(idx_config[idx_idx, "IDX_ID"]), "LATEX-NAME"], " (MB)", sep="", collapse="")
	)
}
    
for(tc_idx in 1:nrow(tc_config)) {
	tc_name = as.character(tc_config[tc_idx, "TC_ID"])

	header = paste(header, c(
		tc_config[as.character(tc_config[tc_idx, "TC_ID"]), "LATEX-NAME"]
	), sep=" & ")

	data<-raw[raw$TC_ID==tc_name,][1:1,]

	table1 = paste(table1, c(
		data[["ALPHABET_SIZE"]],
		sprintf("%.1f", data[["DATA_SIZE"]] / (1024 * 1024)),
		sprintf("%.1f", data[["TREE_SIZE"]] / (1024 * 1024))
	), sep=" & ")

	table2_row = c()
	for(idx_idx in 1:nrow(idx_config)) {
		idx_name = as.character(idx_config[idx_idx, "IDX_ID"])
		data<-raw[raw$TC_ID==tc_name & raw$IDX_ID==idx_name,]
		table2_row = c(table2_row, sprintf("%.2f", data[["CST_SIZE"]] / (1024 * 1024)))
	}
	table2 = paste(table2, table2_row, sep=" & ")
}

cat(header, "\\hline", sep="\\\\\n")
cat(table1, "\\hline", sep="\\\\\n")
cat(table2, "\\hline", sep="\\\\\n")
cat("\\end{tabular}\n")

sink("fig-time.tex")

cat("\\begin{tabular}{ll|*{", nrow(tc_config), "}{r}}\n", sep="")
cat("\\hline\n")

header = c("Operation & CST")
for(tc_idx in 1:nrow(tc_config)) {
    header = paste(header, c(
		tc_config[as.character(tc_config[tc_idx, "TC_ID"]), "LATEX-NAME"]
	), sep=" & ")
}
cat(header, "\\hline", sep="\\\\\n")

format_string = "%.2f"

for(op_idx in 1:nrow(operations)) {
	op_name = as.character(operations[op_idx, "OP_ID"])
	
	for(idx_idx in 1:nrow(idx_config)) {
		idx_name = as.character(idx_config[idx_idx, "IDX_ID"])
		
		if(idx_idx == 1) {
			table = c(paste(as.character(operations[op_idx, "LATEX-NAME"]),
							as.character(idx_config[idx_idx, "LATEX-NAME"]), sep=" & ", collapse=""))
		} else {
			table = c(paste("", as.character(idx_config[idx_idx, "LATEX-NAME"]), sep=" & ", collapse=""))
		}
		
		for(tc_idx in 1:nrow(tc_config)) {
			tc_name = as.character(tc_config[tc_idx, "TC_ID"])

			data<-raw[raw$IDX_ID==idx_name & raw$TC_ID==tc_name,]
			table = paste(table, c(
				sprintf(format_string, data[[paste("CST_", op_name, "_U_TIME", sep="", collapse="")]] / 1000)
			), sep=" & ")
		}

		if(idx_idx == nrow(idx_config)) {
			cat(table, "\\\\\\hline\n", sep="")
		} else {
			cat(table, "\\\\\n", sep="")
		}
	}
}

cat("\\end{tabular}\n")

sink(NULL)
