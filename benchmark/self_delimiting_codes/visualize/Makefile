# Makefile will also automatically generate a tex file with system information

#utility
empty:=
space:= $(empty) $(empty)

auto:	self_delimiting_codes.pdf

self_delimiting_codes.pdf:	../results/result.csv ../results/tc.csv ../results/vat.csv
	$(eval CPUINFO := $(strip $(patsubst \\%, $(space),\
		$(shell cat /proc/cpuinfo | grep "model name.*" | uniq | cut -d':' -f 2))))
	$(eval MEMINFO := $(strip \
		$(firstword $(shell free -k | grep "Mem.*" | uniq | cut -d':' -f 2))))
	$(eval MEMINFO := $(addsuffix " KB", $(MEMINFO)))
	$(eval DISTINFO := $(strip $(patsubst \\%, $(space),\
		$(shell cat /etc/issue))))
	@echo "Creating system information file"
	@echo "\\\\begin{tabular}{ll}" > sysinfo.tex
	@echo "\\\\toprule" >> sysinfo.tex
	@echo "CPU          & \\\\verb\\\\$(CPUINFO)\\\\ \\\\\\\\" >> sysinfo.tex
	@echo "Total Memory & \\\\verb\\\\$(MEMINFO)\\\\ \\\\\\\\" >> sysinfo.tex
	@echo "Distribution & \\\\verb\\\\$(DISTINFO)\\\\ \\\\\\\\" >> sysinfo.tex
	@echo "\\\\bottomrule" >> sysinfo.tex
	@echo "\\\\end{tabular}" >> sysinfo.tex
	@echo "Use pdflatex to generate self_delimiting_codes.pdf"
	@pdflatex self_delimiting_codes.tex >> LaTeX.log 2>&1

cleanall:
	rm -f self_delimiting_codes.pdf self_delimiting_codes.aux \
	self_delimiting_codes.log LaTeX.log sysinfo.tex
	
