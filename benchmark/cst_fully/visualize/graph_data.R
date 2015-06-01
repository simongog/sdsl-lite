source("../../basic_functions.R")

tc_config <- readConfig("../test_case.config", c("TC_ID","PATH","LATEX-NAME","URL"))

operations = data.frame(id=c("LCA", "LETTER", "SLINK", "CHILD", "DEPTH", "PARENT"),
						caption=c("LCA", "Letter", "SLink", "Child", "SDepth", "Parent"),
						filename=c("lca", "letter", "slink", "child", "depth", "parent"))
rownames(operations)<-operations[["OP_ID"]]
colnames(operations)<-c("OP_ID", "LATEX-NAME", "FILE-NAME")

idx_types = data.frame(id=c("FULLY", "FULLY_SDS", "FULLY_BLIND"),
					   caption=c("cst\\_fully", "cst\\_fully\\_sds", "cst\\_fully\\_blind"))
rownames(idx_types)<-idx_types[["IDX_ID"]]
colnames(idx_types)<-c("IDX_ID", "LATEX-NAME")

sampling_rates = c(4, 8, 16, 32, 64, 128)

raw <- data_frame_from_key_value_pairs("../results/all.txt")

for(op_idx in 1:nrow(operations)) {
	op_name = as.character(operations[op_idx, "OP_ID"])
	
	sink(paste("output/", as.character(operations[op_idx, "FILE-NAME"]), ".dat", sep="", collapse=""))
	
	for(delta in sampling_rates) {
		cat(sprintf("%d", delta))
		for(idx_idx in 1:nrow(idx_types)) {
			idx_name = paste(as.character(idx_types[idx_idx, "IDX_ID"]), delta, sep="_")
			
			data<-raw[raw$IDX_ID==idx_name,]
			cat("", as.character(data[["FCST_SIZE"]]))
			cat("", sprintf("%d", data[[paste("FCST_", op_name, "_U_TIME", sep="", collapse="")]]))
		}
		cat("\n")
	}
}

sink(NULL)

