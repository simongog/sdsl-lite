include ../../Make.helper
CFLAGS = $(MY_CXX_FLAGS)
SRC_DIR = src
BIN_DIR = bin
LIBS = -lsdsl

C_OPTIONS:=$(call config_ids,compile_options.config)
TC_IDS:=$(call config_ids,test_case.config)
K2_IDS:=$(call config_ids,k2tree.config)

DL = $(foreach TC_ID,$(TC_IDS),\
		$(call config_select,test_case.config,$(TC_ID),2))
DL_K2T = $(foreach TC_ID,$(TC_IDS),\
		 	$(foreach K2_ID,$(K2_IDS),\
		 		../tmp/K2T.$(TC_ID).$(K2_ID)))

K2_EXECS = $(foreach K2_ID,$(K2_IDS),$(BIN_DIR)/build_$(K2_ID))

RES_FILES = $(foreach K2_ID,$(K2_IDS),\
				$(foreach TC_ID,$(TC_IDS),\
					results/$(K2_ID).$(TC_ID)))

GEN_FILES = $(foreach K2_ID,$(K2_IDS),\
					$(BIN_DIR)/gen_k2_$(K2_ID))

K2T_RES_FILES = $(foreach K2_ID,$(K2_IDS),\
					$(foreach TC_ID,$(TC_IDS),\
						results-k2t/$(K2_ID).$(TC_ID)))

RESULT_FILE=results/all.txt

all: execs

execs: $(K2_EXECS)

timing: execs $(RES_FILES)
	@cat $(RES_FILES) > $(RESULT_FILE)
	@cd visualize;make

../tmp/K2T.%: $(DL) $(GEN_FILES)
	$(eval TC_ID:=$(call dim,1,$*))
	$(eval K2_ID:=$(call dim,2,$*))
	$(eval TC_PATH:=$(call config_select,test_case.config,$(TC_ID),2))
	$(eval NUM_BYTE:=$(call config_select,test_case.config,$(TC_ID),5))
	@$(BIN_DIR)/gen_k2_$(K2_ID) $(TC_PATH) "../tmp/K2T.$(TC_ID).VECTOR" "../tmp/K2T.$(TC_ID).$(K2_ID)"

# Execute $(BIN_DIR)/build_[K2_ID] and write result
results/%: test_case.config $(DL) $(DL_K2T) execs
	$(eval K2_ID:=$(call dim,1,$*))
	$(eval TC_ID:=$(call dim,2,$*))
	$(eval K2_TEX_NAME:=$(call config_select,k2tree.config,$(K2_ID),3))
	$(eval TC_TEX_NAME:=$(call config_select,test_case.config,$(TC_ID),3))
	$(eval K2T:=$(strip $(call config_select,test_case.config,$(TC_ID),6)))
	$(eval TC_PATH:=$(call config_select,test_case.config,$(TC_ID),2))
	$(eval TC_TYPE:=$(call config_select,test_case.config,$(TC_ID),5))
	$(eval TC_SIZE:=$(shell wc -c <$(TC_PATH)))
	$(eval ARGS:="../tmp/K2T.$(TC_ID).$(K2_ID)")
	@echo "Running bin/build_$(K2_ID) on $(TC_ID)"
	@echo "# K2_ID = $(K2_ID)" > $@
	@echo "# TC_ID = $(TC_ID)" >> $@
	@echo "# K2_TEX_NAME = $(K2_TEX_NAME)">>$@
	@echo "# TC_TEX_NAME = $(TC_TEX_NAME)">>$@
	@echo "# TC_SIZE = $(TC_SIZE)">>$@
	@$(BIN_DIR)/build_$(K2_ID) $(ARGS) >> $@

# $(BIN_DIR)/build_[K2_ID]
$(BIN_DIR)/build_%: $(SRC_DIR)/k2_time_and_space.cpp k2tree.config
	$(eval K2_ID:=$(call dim,1,$*))
	$(eval K2_TYPE:=$(call config_select,k2tree.config,$(K2_ID),2))
	@$(MY_CXX) $(CFLAGS) $(C_OPTIONS) -DK2_TYPE="$(K2_TYPE)" -L$(LIB_DIR)\
		$(SRC_DIR)/k2_time_and_space.cpp -I$(INC_DIR) -o $@ $(LIBS)

$(BIN_DIR)/gen_k2_%: $(SRC_DIR)/gen_k2t.cpp
	$(eval K2_ID:=$(call dim,1,$*))
	$(eval K2_TYPE:=$(call config_select,k2tree.config,$(K2_ID),2))
	@$(MY_CXX) $(CFLAGS) $(C_OPTIONS) -DK2_TYPE="$(K2_TYPE)" -L$(LIB_DIR)\
		$(SRC_DIR)/gen_k2t.cpp -I$(INC_DIR) -o $@ $(LIBS) -ldivsufsort -ldivsufsort64

include ../Make.download

clean-build:
	@echo "Remove executables"
	rm -rf $(BIN_DIR)/build*
	rm -rf $(BIN_DIR)/gen*

clean-result:
	@echo "Remove results"
	rm -rf results/*

cleanall: clean-build clean-result
