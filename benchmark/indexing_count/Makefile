include ../Make.helper
CFLAGS = $(MY_CXX_FLAGS) # in compile_options.config
LIBS = -lsdsl -ldivsufsort -ldivsufsort64
SRC_DIR = src
TMP_DIR = ../tmp
PAT_DIR = pattern
BIN_DIR = bin

TC_PATHS:=$(call config_column,test_case.config,2)
TC_IDS:=$(call config_ids,test_case.config)
IDX_IDS:=$(call config_ids,index.config)
COMPILE_IDS:=$(call config_ids,compile_options.config)

RESULT_FILE=results/all.txt

QUERY_EXECS = $(foreach IDX_ID,$(IDX_IDS),\
		         $(foreach COMPILE_ID,$(COMPILE_IDS),$(BIN_DIR)/query_idx_$(IDX_ID).$(COMPILE_ID)))
BUILD_EXECS = $(foreach IDX_ID,$(IDX_IDS),$(BIN_DIR)/build_idx_$(IDX_ID))
INFO_EXECS  = $(foreach IDX_ID,$(IDX_IDS),$(BIN_DIR)/info_$(IDX_ID))
PATTERNS    = $(foreach TC_ID,$(TC_IDS),$(PAT_DIR)/$(TC_ID).pattern)
INDEXES     = $(foreach IDX_ID,$(IDX_IDS),\
		         $(foreach TC_ID,$(TC_IDS),indexes/$(TC_ID).$(IDX_ID)))
INFO_FILES  = $(foreach IDX_ID,$(IDX_IDS),\
		         $(foreach TC_ID,$(TC_IDS),info/$(TC_ID).$(IDX_ID).html))
TIME_FILES  = $(foreach IDX_ID,$(IDX_IDS),\
		         $(foreach TC_ID,$(TC_IDS),\
					 $(foreach COMPILE_ID,$(COMPILE_IDS),results/$(TC_ID).$(IDX_ID).$(COMPILE_ID))))
COMP_FILES  = $(addsuffix .z.info,$(TC_PATHS))

all: $(BUILD_EXECS) $(QUERY_EXECS) $(INFO_EXECS)

info: $(INFO_EXECS) $(INFO_FILES)

indexes: $(INDEXES)

input: $(TC_PATHS)

pattern: input $(PATTERNS) $(BIN_DIR)/genpatterns

compression: input $(COMP_FILES)

timing: input $(INDEXES) pattern $(TIME_FILES) compression info
	@cat $(TIME_FILES) > $(RESULT_FILE)
	@cd visualize; make

# results/[TC_ID].[IDX_ID].[COMPILE_ID]
results/%: $(BUILD_EXECS) $(QUERY_EXECS) $(PATTERNS) $(INDEXES)
	$(eval TC_ID:=$(call dim,1,$*))
	$(eval IDX_ID:=$(call dim,2,$*))
	$(eval COMPILE_ID:=$(call dim,3,$*))
	$(eval TC_NAME:=$(call config_select,test_case.config,$(TC_ID),3))
	@echo "# TC_ID = $(TC_ID)" >> $@
	@echo "# IDX_ID = $(IDX_ID)" >> $@
	@echo "# COMPILE_ID = $(COMPILE_ID)" >> $@
	@echo "# test_case = $(TC_NAME)" >>  $@
	@echo "Run timing for $(IDX_ID).$(COMPILE_ID) on $(TC_ID)"
	$(BIN_DIR)/query_idx_$(IDX_ID).$(COMPILE_ID) \
		indexes/$(TC_ID) C < $(PAT_DIR)/$(TC_ID).pattern 2>> $@


# indexes/[TC_ID].[IDX_ID]
indexes/%: $(BUILD_EXECS)
	$(eval TC_ID:=$(call dim,1,$*))
	$(eval IDX_ID:=$(call dim,2,$*))
	$(eval TC:=$(call config_select,test_case.config,$(TC_ID),2))
	@echo "Building index $(IDX_ID) on $(TC)"
	$(BIN_DIR)/build_idx_$(IDX_ID) $(TC) $(TMP_DIR) $@

# info/[TC_ID].[IDX_ID]
info/%.html: $(INDEXES)
	$(eval TC_ID:=$(call dim,1,$*))
	$(eval IDX_ID:=$(call dim,2,$*))
	@echo "Generating info for $(IDX_ID) on $(TC_ID)"
	@$(BIN_DIR)/info_$(IDX_ID) indexes/$(TC_ID).$(IDX_ID) > $@

$(PAT_DIR)/%.pattern: $(BIN_DIR)/genpatterns
	@echo $*
	$(eval TC:=$(call config_select,test_case.config,$*,2))
	$(BIN_DIR)/genpatterns $(TC) 20 50000 $@

$(BIN_DIR)/genpatterns: $(SRC_DIR)/genpatterns.c
	@echo "Build pattern generation program"
	@$(MY_CC) -O3 -o $@ $(SRC_DIR)/genpatterns.c

# $(BIN_DIR)/build_idx_[IDX_ID]
$(BIN_DIR)/build_idx_%: $(SRC_DIR)/build_index_sdsl.cpp index.config
	$(eval IDX_TYPE:=$(call config_select,index.config,$*,2))
	@echo "Compiling build_idx_$*"
	@$(MY_CXX) $(CFLAGS) -O3 -DNDEBUG \
					-DSUF=\"$*\" -DCSA_TYPE="$(IDX_TYPE)" \
					-L$(LIB_DIR) $(SRC_DIR)/build_index_sdsl.cpp \
					-I$(INC_DIR) -o $@ $(LIBS)

# Targets for the count experiment. $(BIN_DIR)/count_queries_[IDX_ID].[COMPILE_ID]
$(BIN_DIR)/query_idx_%: $(SRC_DIR)/run_queries_sdsl.cpp index.config
	$(eval IDX_ID:=$(call dim,1,$*))
	$(eval COMPILE_ID:=$(call dim,2,$*))
	$(eval IDX_TYPE:=$(call config_select,index.config,$(IDX_ID),2))
	$(eval COMPILE_OPTIONS:=$(call config_select,compile_options.config,$(COMPILE_ID),2))
	@echo "Compiling query_idx_$*"
	@$(MY_CXX) $(CFLAGS) $(COMPILE_OPTIONS) \
					-DSUF=\"$(IDX_ID)\" -DCSA_TYPE="$(IDX_TYPE)" \
					-L$(LIB_DIR) $(SRC_DIR)/run_queries_sdsl.cpp \
					-I$(INC_DIR) -o $@ $(LIBS)

# Targets for the executables which output the indexes structure.
$(BIN_DIR)/info_%: $(SRC_DIR)/info.cpp index.config
	$(eval IDX_TYPE:=$(call config_select,index.config,$*,2))
	@echo "Compiling info_$*"
	@$(MY_CXX) $(CFLAGS) -O3 -DNDEBUG \
					-DSUF=\"$*\" -DCSA_TYPE="$(IDX_TYPE)" \
					-L$(LIB_DIR) $(SRC_DIR)/info.cpp \
					-I$(INC_DIR) -o $@ $(LIBS)

include ../Make.download

clean-build:
	@echo "Remove executables"
	@rm -f $(QUERY_EXECS) $(BUILD_EXECS) $(INFO_EXECS)

clean:
	@echo "Remove executables and indexes"
	@rm -f $(QUERY_EXECS) $(BUILD_EXECS) $(INFO_EXECS) \
		   $(INFO_FILES) $(INDEXES) $(BIN_DIR)/genpatterns

cleanresults:
	@echo "Remove result files"
	@rm -f $(TIME_FILES) $(RESULT_FILE) $(INFO_FILES)
	@rm -f $(PATTERNS)

cleanall: clean cleanresults
	@echo "Remove all generated files."
	@rm -f $(TMP_DIR)/*
	@rm -f $(PAT_DIR)/*
