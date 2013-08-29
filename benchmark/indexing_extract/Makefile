include ../Make.helper
CFLAGS = $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) 
LIBS = -lsdsl -ldivsufsort -ldivsufsort64
SRC_DIR = src
TMP_DIR = ../tmp
IVL_DIR = intervals
BIN_DIR = bin

TC_PATHS:=$(call config_column,test_case.config,2)
TC_IDS:=$(call config_ids,test_case.config)
IDX_IDS:=$(call config_ids,index.config)
SAMPLE_IDS:=$(call config_ids,sample.config)

RESULT_FILE=results/all.txt

QUERY_EXECS = $(foreach IDX_ID,$(IDX_IDS),\
		        $(foreach SAMPLE_ID,$(SAMPLE_IDS),$(BIN_DIR)/query_idx_$(IDX_ID).$(SAMPLE_ID)))
BUILD_EXECS = $(foreach IDX_ID,$(IDX_IDS),\
				$(foreach SAMPLE_ID,$(SAMPLE_IDS),$(BIN_DIR)/build_idx_$(IDX_ID).$(SAMPLE_ID)))
INFO_EXECS  = $(foreach IDX_ID,$(IDX_IDS),\
				$(foreach SAMPLE_ID,$(SAMPLE_IDS),$(BIN_DIR)/info_$(IDX_ID).$(SAMPLE_ID)))
INTERVALS   = $(foreach TC_ID,$(TC_IDS),$(IVL_DIR)/$(TC_ID).interval)
INDEXES     = $(foreach IDX_ID,$(IDX_IDS),\
		        $(foreach TC_ID,$(TC_IDS),\
				  $(foreach SAMPLE_ID,$(SAMPLE_IDS),indexes/$(TC_ID).$(IDX_ID).$(SAMPLE_ID))))
INFO_FILES  = $(foreach IDX_ID,$(IDX_IDS),\
		      	$(foreach TC_ID,$(TC_IDS),\
				  $(foreach SAMPLE_ID,$(SAMPLE_IDS),info/$(TC_ID).$(IDX_ID).$(SAMPLE_ID).json)))
TIME_FILES  = $(foreach IDX_ID,$(IDX_IDS),\
		        $(foreach TC_ID,$(TC_IDS),\
				  $(foreach SAMPLE_ID,$(SAMPLE_IDS),results/$(TC_ID).$(IDX_ID).$(SAMPLE_ID))))
COMP_FILES  = $(addsuffix .z.info,$(TC_PATHS))

all: $(BUILD_EXECS) $(QUERY_EXECS) $(INFO_EXECS)

info: $(INFO_EXECS) $(INFO_FILES) 

indexes: $(INDEXES)

intervals: input $(INTERVALS) $(BIN_DIR)/genintervals

input: $(TC_PATHS)

compression: input $(COMP_FILES)

timing: input $(INDEXES) intervals $(TIME_FILES) compression info
	@cat $(TIME_FILES) > $(RESULT_FILE)
	@cd visualize; make

# results/[TC_ID].[IDX_ID].[SAMPLE_ID]
results/%: $(BUILD_EXECS) $(QUERY_EXECS) $(PATTERNS) $(INDEXES)
	$(eval TC_ID:=$(call dim,1,$*)) 
	$(eval IDX_ID:=$(call dim,2,$*)) 
	$(eval SAMPLE_ID:=$(call dim,3,$*)) 
	$(eval TC_NAME:=$(call config_select,test_case.config,$(TC_ID),3))
	$(eval S_SA:=$(call config_select,sample.config,$(SAMPLE_ID),2))
	$(eval S_ISA:=$(call config_select,sample.config,$(SAMPLE_ID),3))
	@echo "# TC_ID = $(TC_ID)" >> $@
	@echo "# IDX_ID = $(IDX_ID)" >> $@
	@echo "# test_case = $(TC_NAME)" >>  $@
	@echo "# SAMPLE_ID = $(SAMPLE_ID)" >> $@
	@echo "# S_SA = $(S_SA)" >> $@
	@echo "# S_ISA = $(S_ISA)" >> $@
	@echo "Run timing for $(IDX_ID).$(SAMPLE_ID) on $(TC_ID)"
	@$(BIN_DIR)/query_idx_$(IDX_ID).$(SAMPLE_ID) \
		indexes/$(TC_ID) E < $(IVL_DIR)/$(TC_ID).interval 2>> $@ 
 

# indexes/[TC_ID].[IDX_ID].[SAMPLE_ID]
indexes/%: $(BUILD_EXECS)
	$(eval TC_ID:=$(call dim,1,$*)) 
	$(eval IDX_ID:=$(call dim,2,$*)) 
	$(eval SAMPLE_ID:=$(call dim,3,$*))
	$(eval TC:=$(call config_select,test_case.config,$(TC_ID),2))
	@echo "Building index $(IDX_ID).$(SAMPLE_ID) on $(TC)"
	@$(BIN_DIR)/build_idx_$(IDX_ID).$(SAMPLE_ID) $(TC) $(TMP_DIR) $@

# info/[TC_ID].[IDX_ID].[SAMPLE_ID]
info/%.json: $(INDEXES)
	$(eval TC_ID:=$(call dim,1,$*)) 
	$(eval IDX_ID:=$(call dim,2,$*)) 
	$(eval SAMPLE_ID:=$(call dim,3,$*)) 
	@echo "Generating info for $(IDX_ID) on $(TC_ID)"
	$(BIN_DIR)/info_$(IDX_ID).$(SAMPLE_ID) indexes/$(TC_ID).$(IDX_ID).$(SAMPLE_ID) > $@ 

# $(IVL_DIR)/[TC_ID].interval
$(IVL_DIR)/%.interval: $(BIN_DIR)/genintervals
	@echo "Generating intervals for $*"
	$(eval TC:=$(call config_select,test_case.config,$*,2))
	@$(BIN_DIR)/genintervals $(TC) 512 10000 $@ 2> /dev/null

$(BIN_DIR)/genintervals: $(SRC_DIR)/genintervals.c
	@echo "Build interval generation program"
	@$(MY_CC) -O3 -o $@ $(SRC_DIR)/genintervals.c 


# $(BIN_DIR)/build_idx_[IDX_ID].[SAMPLE_ID]
$(BIN_DIR)/build_idx_%: $(SRC_DIR)/build_index_sdsl.cpp index.config sample.config
	$(eval IDX_ID:=$(call dim,1,$*)) 
	$(eval SAMPLE_ID:=$(call dim,2,$*)) 
	$(eval IDX_TYPE:=$(call config_select,index.config,$(IDX_ID),2))
	$(eval S_SA:=$(call config_select,sample.config,$(SAMPLE_ID),2))
	$(eval S_ISA:=$(call config_select,sample.config,$(SAMPLE_ID),3))
	$(eval IDX_TYPE:=$(subst S_SA,$(S_SA),$(IDX_TYPE)))
	$(eval IDX_TYPE:=$(subst S_ISA,$(S_ISA),$(IDX_TYPE)))
	@echo "Compiling build_idx_$*"
	@$(MY_CXX) $(CFLAGS) -DSUF=\"$*\" -DCSA_TYPE="$(IDX_TYPE)" \
					-DS_SA=$(S_SA) -DS_ISA=$(S_ISA) \
					-L$(LIB_DIR) $(SRC_DIR)/build_index_sdsl.cpp \
					-I$(INC_DIR) -o $@ $(LIBS)

# Targets for the count experiment. $(BIN_DIR)/count_queries_[IDX_ID].[SAMPLE_ID]
$(BIN_DIR)/query_idx_%: $(SRC_DIR)/run_queries_sdsl.cpp index.config sample.config
	$(eval IDX_ID:=$(call dim,1,$*)) 
	$(eval SAMPLE_ID:=$(call dim,2,$*)) 
	$(eval IDX_TYPE:=$(call config_select,index.config,$(IDX_ID),2))
	$(eval S_SA:=$(call config_select,sample.config,$(SAMPLE_ID),2))
	$(eval S_ISA:=$(call config_select,sample.config,$(SAMPLE_ID),3))
	$(eval IDX_TYPE:=$(subst S_SA,$(S_SA),$(IDX_TYPE)))
	$(eval IDX_TYPE:=$(subst S_ISA,$(S_ISA),$(IDX_TYPE)))
	@echo "Compiling query_idx_$*"
	@$(MY_CXX) $(CFLAGS) -DSUF="$(IDX_ID).$(SAMPLE_ID)" -DCSA_TYPE="$(IDX_TYPE)" \
			         -L$(LIB_DIR) $(SRC_DIR)/run_queries_sdsl.cpp \
			         -I$(INC_DIR) -o $@ $(LIBS)

# Targets for the info executable. $(BIN_DIR)/info_[IDX_ID].[SAMPLE_ID]
$(BIN_DIR)/info_%: $(SRC_DIR)/info.cpp index.config 
	$(eval IDX_ID:=$(call dim,1,$*)) 
	$(eval SAMPLE_ID:=$(call dim,2,$*)) 
	$(eval IDX_TYPE:=$(call config_select,index.config,$(IDX_ID),2))
	$(eval S_SA:=$(call config_select,sample.config,$(SAMPLE_ID),2))
	$(eval S_ISA:=$(call config_select,sample.config,$(SAMPLE_ID),3))
	$(eval IDX_TYPE:=$(subst S_SA,$(S_SA),$(IDX_TYPE)))
	$(eval IDX_TYPE:=$(subst S_ISA,$(S_ISA),$(IDX_TYPE)))
	@echo "Compiling info_$*"
	@$(MY_CXX) $(CFLAGS) -DSUF=\"$*\" -DCSA_TYPE="$(IDX_TYPE)" \
					-L$(LIB_DIR) $(SRC_DIR)/info.cpp \
					-I$(INC_DIR) -o $@ $(LIBS)

include ../Make.download

clean-build:
	@echo "Remove executables"
	@rm -f $(QUERY_EXECS) $(BUILD_EXECS) $(INFO_EXECS)

clean:
	@echo "Remove executables"
	@rm -f $(QUERY_EXECS) $(BUILD_EXECS) $(INFO_EXECS) \
		   $(BIN_DIR)/genintervals

cleanresults: 
	@echo "Remove result files"
	@rm -f $(TIME_FILES) $(RESULT_FILE)

cleanall: clean cleanresults
	@echo "Remove all generated files."
	@rm -f $(INDEXES) $(INFO_FILES) $(PATTERNS)
	@rm -f $(TMP_DIR)/* 
	@rm -f $(INTERVALS)
