include ../../Make.helper
CFLAGS = $(MY_CXX_FLAGS) 
SRC_DIR = src
BIN_DIR = bin
LIBS = -lsdsl 

C_OPTIONS:=$(call config_ids,compile_options.config)
TC_IDS:=$(call config_ids,test_case.config)
WT_IDS:=$(call config_ids,wt.config)

DL = $(foreach TC_ID,$(TC_IDS),\
		$(call config_select,test_case.config,$(TC_ID),2))
DL_BWT = $(foreach TC_ID,$(TC_IDS),../tmp/BWT.$(TC_ID))

WT_EXECS = $(foreach WT_ID,$(WT_IDS),$(BIN_DIR)/build_$(WT_ID))

RES_FILES = $(foreach WT_ID,$(WT_IDS),\
				$(foreach TC_ID,$(TC_IDS),\
					results/$(WT_ID).$(TC_ID)))

BWT_RES_FILES = $(foreach WT_ID,$(WT_IDS),\
					$(foreach TC_ID,$(TC_IDS),\
						results-bwt/$(WT_ID).$(TC_ID)))

RESULT_FILE=results/all.txt

all: execs 

execs: $(WT_EXECS)

timing: execs $(RES_FILES) 
	@cat $(RES_FILES) > $(RESULT_FILE)
	@cd visualize;make

../tmp/BWT.%: $(DL) $(BIN_DIR)/gen_bwt 
	$(eval TC_ID:=$*)
	$(eval TC_PATH:=$(call config_select,test_case.config,$(TC_ID),2))
	$(eval BWT:=$(strip $(call config_select,test_case.config,$(TC_ID),6)))
	$(eval NUM_BYTE:=$(call config_select,test_case.config,$(TC_ID),5))
	@$(BIN_DIR)/gen_bwt $(TC_PATH) ../tmp/BWT.$(TC_ID) ../tmp $(BWT) $(NUM_BYTE)

# Execute $(BIN_DIR)/build_[WT_ID] and write result
results/%: test_case.config $(DL) $(DL_BWT) execs 
	$(eval WT_ID:=$(call dim,1,$*))
	$(eval TC_ID:=$(call dim,2,$*))
	$(eval WT_TEX_NAME:=$(call config_select,wt.config,$(WT_ID),3))
	$(eval TC_TEX_NAME:=$(call config_select,test_case.config,$(TC_ID),3))
	$(eval BWT:=$(strip $(call config_select,test_case.config,$(TC_ID),6)))
	$(eval TC_PATH:=$(if $(findstring BWT,$(BWT)),"../tmp/BWT.$(TC_ID)",$(call config_select,test_case.config,$(TC_ID),2)))
	$(eval TC_TYPE:=$(call config_select,test_case.config,$(TC_ID),5))
	$(eval TC_SIZE:=$(shell wc -c <$(TC_PATH)))
	$(eval ARGS:=$(TC_PATH) $(TC_TYPE) $(TC_ID) $(WT_ID))
	@echo "Running bin/build_$(WT_ID) on $(TC_ID)"
	@echo "# WT_ID = $(WT_ID)" > $@
	@echo "# TC_ID = $(TC_ID)" >> $@
	@echo "# WT_TEX_NAME = $(WT_TEX_NAME)">>$@
	@echo "# TC_TEX_NAME = $(TC_TEX_NAME)">>$@
	@echo "# TC_SIZE = $(TC_SIZE)">>$@
	@$(BIN_DIR)/build_$(WT_ID) $(ARGS) >> $@	

# $(BIN_DIR)/build_[WT_ID]
$(BIN_DIR)/build_%: $(SRC_DIR)/wt_time_and_space.cpp wt.config 
	$(eval WT_ID:=$(call dim,1,$*))
	$(eval WT_TYPE:=$(call config_select,wt.config,$(WT_ID),2))
	@echo "Compiling build_$*"
	@$(MY_CXX) $(CFLAGS) $(C_OPTIONS) -DWT_TYPE="$(WT_TYPE)" -L$(LIB_DIR)\
		$(SRC_DIR)/wt_time_and_space.cpp -I$(INC_DIR) -o $@ $(LIBS)

$(BIN_DIR)/gen_bwt: $(SRC_DIR)/gen_bwt.cpp
	@$(MY_CXX) $(CFLAGS) $(C_OPTIONS) -L$(LIB_DIR)\
		$(SRC_DIR)/gen_bwt.cpp -I$(INC_DIR) -o $@ $(LIBS) -ldivsufsort -ldivsufsort64

include ../Make.download

clean-build:
	@echo "Remove executables"
	rm -f $(BIN_DIR)/build*

clean-result:
	@echo "Remove results"
	rm -f results/*

cleanall: clean-build clean-result 
