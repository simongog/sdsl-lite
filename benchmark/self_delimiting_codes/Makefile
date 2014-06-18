include ../../Make.helper
SRC_DIR = src
BIN_DIR = bin
LIBS = -lsdsl
RES_FILE = results/result.csv	#result file of benchmark
VAT_FILE = results/vat.csv	#vector assignment table (vector name -> sdsl type)
TC_FILE = results/tc.csv	#test case table (contains only test case names)

#utility
empty:=
space:= $(empty) $(empty)
comma:= ,

#load test cases
TC_IDS   := $(call config_ids,test_case.config)
TC_FILES := $(foreach TC_ID,$(TC_IDS),\
		$(call config_select,test_case.config,$(TC_ID),2))

all:	$(RES_FILE)

timing: $(RES_FILE)
	@cd visualize;make


#compilation and creation of vector assignment table
$(BIN_DIR)/sdcbenchmark: $(SRC_DIR)/sdc_benchmark.cpp vectors.config compile_options.config
	$(eval VTYPES := $(subst $(space),$(comma),$(strip $(call config_column,vectors.config,2))))
	$(eval VNAMES := $(subst $(space),\"$(comma)\",$(strip $(call config_column,vectors.config,3))))
	$(eval VNAMES := $(addprefix {\",$(VNAMES)))
	$(eval VNAMES := $(addsuffix \"},$(VNAMES)))
	$(eval C_OPTIONS:=$(call config_ids,compile_options.config))
	@echo "Compiling build for vectors $(VNAMES)"
	@$(MY_CXX) $(MY_CXX_FLAGS) $(C_OPTIONS) -DVTYPES="$(VTYPES)" -DVNAMES="$(VNAMES)" -L$(LIB_DIR)\
		"$(SRC_DIR)/sdc_benchmark.cpp" -I$(INC_DIR) -o "$(BIN_DIR)/sdcbenchmark" $(LIBS)
	$(eval V_IDS := $(call config_ids,vectors.config))
	$(eval V_ASSIGNMENTTABLE := $(subst $(space),\n,$(strip $(foreach V_ID,$(V_IDS),\
			$(call config_select,vectors.config,$(V_ID),3);$(call config_select,vectors.config,$(V_ID),2)))))
	@echo "Writing Vector Assignment Table"
	@echo "vector;sdsltype" > $(VAT_FILE)
	@echo "$(V_ASSIGNMENTTABLE)" >> $(VAT_FILE)
	

#execution and creation of test case table
$(RES_FILE):	test_case.config $(TC_FILES) $(BIN_DIR)/sdcbenchmark
	$(eval ARGS := $(foreach TC_ID,$(TC_IDS),\
		$(call config_select,test_case.config,$(TC_ID),3) $(space) \
		$(call config_select,test_case.config,$(TC_ID),2) $(space) \
		$(call config_select,test_case.config,$(TC_ID),5) ) )
	@echo "Executing Benchmark"
	@$(BIN_DIR)/sdcbenchmark $(ARGS) | tee $(RES_FILE)
	$(eval TC_TABLE := $(subst $(space),\n,$(strip $(call config_column,test_case.config,3))))
	@echo "Writing Test Case file"
	@echo "testcase\\nOverall" > $(TC_FILE)	
	@echo "$(TC_TABLE)" >> $(TC_FILE)		

include ../Make.download

clean-build:
	@echo "Remove executables"
	rm -f $(BIN_DIR)/sdcbenchmark

clean-result:
	@echo "Remove results"
	rm -f results/*

cleanall: clean-build clean-result
	@cd visualize;make cleanall
