include ../../Make.helper
CFLAGS = $(MY_CXX_FLAGS) 
SRC_DIR = src
BIN_DIR = bin
LIBS = -lsdsl -ldivsufsort -ldivsufsort64

C_OPTIONS:=$(call config_ids,compile_options.config)
TC_IDS:=$(call config_ids,test_case.config)
LCP_IDS:=$(call config_ids,lcp.config)


DL = ${foreach TC_ID,$(TC_IDS),$(call config_select,test_case.config,$(TC_ID),2)}

LCP_EXECS = $(foreach LCP_ID,$(LCP_IDS),$(BIN_DIR)/build_$(LCP_ID))

RES_FILES = $(foreach TC_ID,$(TC_IDS),\
					results/$(TC_ID))

RESULT_FILE=results/all.txt

execs: $(BIN_DIR)/prep_sa_bwt $(LCP_EXECS)

timing: execs $(RES_FILES)
	@cat $(RES_FILES) > $(RESULT_FILE)
	@cd visualize;make

$(BIN_DIR)/prep_sa_bwt: $(SRC_DIR)/create_sa_bwt.cpp 
	@echo "Compiling prep_sa_bwt"
	@$(MY_CXX) $(CFLAGS) $(C_OPTIONS) -L${LIB_DIR}\
		$(SRC_DIR)/create_sa_bwt.cpp -I${INC_DIR} -o bin/prep_sa_bwt $(LIBS)
	
precalc%: test_case.config $(DL) lcp.config
	$(eval TC_ID:=$(call dim,1,$*))
	$(eval LCP_TEX_NAME:=$(call config_select,lcp.config,$(LCP_ID),3))
	$(eval TC_TEX_NAME:=$(call config_select,test_case.config,$(TC_ID),3))
	$(eval TC_PATH:=$(call config_select,test_case.config,$(TC_ID),2))
	$(eval TC_SIZE:=$(shell wc -c <$(TC_PATH)))
	@echo "Running test case: $(TC_ID)"
	@echo "# TC_ID = $(TC_ID)" > results/$(TC_ID)
	@echo "# TC_TEX_NAME = $(TC_TEX_NAME)">> results/$(TC_ID)
	@echo "# TC_SIZE = $(TC_SIZE)">> results/$(TC_ID)
	@$(BIN_DIR)/prep_sa_bwt $(TC_PATH) >> results/$(TC_ID) 

results/%: precalc%
	@$(foreach LCP_EXEC,$(LCP_EXECS),$(shell $(LCP_EXEC) >>$@;rm -f lcp_tmp.sdsl isa_tmp.sdsl)) 
	@rm *.sdsl

$(BIN_DIR)/build_%: $(SRC_DIR)/create_lcp.cpp lcp.config 
	$(eval LCP_ID:=$(call dim,1,$*))
	$(eval LCP_TYPE:=$(call config_select,lcp.config,$(LCP_ID),2))
	@echo "Compiling build_$*"
	@$(MY_CXX) $(CFLAGS) $(C_OPTIONS) -DLCP_TYPE="$(LCP_TYPE)" -DLCPID="$(LCP_ID)" -L${LIB_DIR}\
		$(SRC_DIR)/create_lcp.cpp -I${INC_DIR} -o $@ $(LIBS)


include ../Make.download

clean-build:
	@echo "Remove executables"
	rm -f $(BIN_DIR)/build*
	rm -f $(BIN_DIR)/prep*

clean-result:
	@echo "Remove results"
	rm -f results/*

cleanall: clean-build clean-result 
