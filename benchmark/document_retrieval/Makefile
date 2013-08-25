include ../Make.helper
CXX_FLAGS = $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(INC_DIR) -L$(LIB_DIR)  
LIBS = -lsdsl -ldivsufsort -ldivsufsort64
SRC_DIR = src
TMP_DIR = ../tmp
PAT_DIR = pattern
BIN_DIR = bin

TC_PATHS:=$(call config_column,test_case.config,2)
TC_IDS:=$(call config_ids,test_case.config)
IDX_IDS:=$(call config_ids,index.config)
PAT_LENS:=$(call config_column,pattern_length.config,1)

RESULT_FILE=results/all.txt

QUERY_EXECS = $(foreach IDX_ID,$(IDX_IDS),$(BIN_DIR)/query_idx_$(IDX_ID))
BUILD_EXECS = $(foreach IDX_ID,$(IDX_IDS),$(BIN_DIR)/build_idx_$(IDX_ID))
SIZES_EXECS = $(foreach IDX_ID,$(IDX_IDS),$(BIN_DIR)/size_of_idx_$(IDX_ID))
PATTERNS    = $(foreach TC_ID,$(TC_IDS),\
				$(foreach PAT_LEN,$(PAT_LENS),$(PAT_DIR)/$(TC_ID).$(PAT_LEN).pattern))
INDEXES     = $(foreach IDX_ID,$(IDX_IDS),\
		         $(foreach TC_ID,$(TC_IDS),indexes/$(TC_ID).$(IDX_ID)))
SIZES       = $(foreach IDX_ID,$(IDX_IDS),\
		         $(foreach TC_ID,$(TC_IDS),info/$(TC_ID).$(IDX_ID).size))
INFO_FILES  = $(foreach IDX_ID,$(IDX_IDS),\
		         $(foreach TC_ID,$(TC_IDS),info/$(TC_ID).$(IDX_ID).json))
TIME_FILES  = $(foreach IDX_ID,$(IDX_IDS),\
		         $(foreach TC_ID,$(TC_IDS),\
					 $(foreach PAT_LEN,$(PAT_LENS),results/$(TC_ID).$(IDX_ID).$(PAT_LEN))))
COMP_FILES  = $(addsuffix .z.info,$(TC_PATHS))
	
all: $(BUILD_EXECS) $(QUERY_EXECS) pattern

info: $(SIZES_EXECS) $(SIZES)
	@cat $(SIZES) > info/sizes.txt

indexes: $(INDEXES)

input: $(TC_PATHS)

pattern: input $(PATTERNS) $(BIN_DIR)/gen_pattern

compression: input $(COMP_FILES)

timing: input $(INDEXES) $(TIME_FILES) pattern compression info
	@cat $(TIME_FILES) > $(RESULT_FILE)
	@cd visualize; make

# results/[TC_ID].[IDX_ID].[PAT_LEN]
results/%: $(BUILD_EXECS) $(QUERY_EXECS) $(INDEXES) $(PATTERNS)
	$(eval TC_ID:=$(call dim,1,$*)) 
	$(eval IDX_ID:=$(call dim,2,$*)) 
	$(eval PAT_LEN:=$(call dim,3,$*)) 
	$(eval TC_NAME:=$(call config_select,test_case.config,$(TC_ID),3))
	@echo "# TC_ID = $(TC_ID)" > $@
	@echo "# IDX_ID = $(IDX_ID)" >> $@
	@echo "# test_case = $(TC_NAME)" >>  $@
	@echo "Run timing for $(IDX_ID) on $(TC_ID) with patterns of length $(PAT_LEN)"
	@$(BIN_DIR)/query_idx_$(IDX_ID) indexes/$(TC_ID).$(IDX_ID) \
		     $(PAT_DIR)/$(TC_ID).$(PAT_LEN).pattern >> $@ 
 

# indexes/[TC_ID].[IDX_ID]
indexes/%: $(BUILD_EXECS)
	$(eval TC_ID:=$(call dim,1,$*)) 
	$(eval IDX_ID:=$(call dim,2,$*)) 
	$(eval TC:=$(call config_select,test_case.config,$(TC_ID),2))
	@echo "Building index $(IDX_ID) on $(TC)"
	@$(BIN_DIR)/build_idx_$(IDX_ID) $(TC) $(TMP_DIR) $@

# info/[TC_ID].[IDX_ID]
info/%.size: $(INDEXES)
	$(eval TC_ID:=$(call dim,1,$*)) 
	$(eval IDX_ID:=$(call dim,2,$*)) 
	$(eval TC:=$(call config_select,test_case.config,$(TC_ID),2))
	$(eval SIZE:=$(call file_size,$(TC)))
	@echo "# TC_ID = $(TC_ID)" > $@
	@echo "# IDX_ID = $(IDX_ID)" >> $@
	@echo "# text_size = $(SIZE)" >> $@
	
	@echo "Get size of index for $(IDX_ID) on $(TC_ID)"
	@$(BIN_DIR)/size_of_idx_$(IDX_ID) indexes/$(TC_ID).$(IDX_ID) >> $@
	
# $(PAT_DIR)/[TC_ID].[PAT_LEN].pattern	
$(PAT_DIR)/%.pattern: $(BIN_DIR)/gen_pattern
	@echo $*
	$(eval TC_ID:=$(call dim,1,$*)) 
	$(eval PAT_LEN:=$(call dim,2,$*)) 
	$(eval TC:=$(call config_select,test_case.config,$(TC_ID),2))
	$(BIN_DIR)/gen_pattern $(TC) $(TMP_DIR)/$(TC_ID).pat.csa $(TMP_DIR) $(PAT_LEN) 200 $@

$(BIN_DIR)/gen_pattern: $(SRC_DIR)/gen_pattern.cpp
	@echo "Build pattern generation program"
	$(MY_CXX) $(CXX_FLAGS) $(SRC_DIR)/gen_pattern.cpp \
				-L$(LIB_DIR) -I$(INC_DIR) -o $@ $(LIBS)

$(BIN_DIR)/size_of_idx_%: $(SRC_DIR)/size_of_idx.cpp
	$(eval IDX_TYPE:=$(call config_select,index.config,$*,2))
	@echo "Build size info program"
	$(MY_CXX) $(CXX_FLAGS) $(SRC_DIR)/size_of_idx.cpp \
				-DIDX_TYPE="$(IDX_TYPE)" \
				-L$(LIB_DIR) -I$(INC_DIR) -o $@ $(LIBS)

# $(BIN_DIR)/build_idx_[IDX_ID]
$(BIN_DIR)/build_idx_%: $(SRC_DIR)/build_idx.cpp index.config
	$(eval IDX_TYPE:=$(call config_select,index.config,$*,2))
	@echo "Compiling build_idx_$*"
	$(MY_CXX) $(CXX_FLAGS) \
					-DIDX_SUF=\"$*\" -DIDX_TYPE="$(IDX_TYPE)" \
					-L$(LIB_DIR) $(SRC_DIR)/build_idx.cpp \
					-I$(INC_DIR) -o $@ $(LIBS)

# Targets for the count experiment. $(BIN_DIR)/query_idx_[IDX_ID]
$(BIN_DIR)/query_idx_%: $(SRC_DIR)/query_idx.cpp index.config 
	$(eval IDX_TYPE:=$(call config_select,index.config,$*,2))
	@echo "Compiling query_idx_$*"
	$(MY_CXX) $(CXX_FLAGS) \
					-DIDX_ID=\"$*\" -DIDX_TYPE="$(IDX_TYPE)" \
					-L$(LIB_DIR) $(SRC_DIR)/query_idx.cpp \
					-I$(INC_DIR) -o $@ $(LIBS)

include ../Make.download

clean:
	@echo "Remove executables and indexes"
	@rm -f $(QUERY_EXECS) $(SIZES_EXECS) $(BUILD_EXECS) \
		   $(INFO_FILES) $(INDEXES) $(BIN_DIR)/gen_pattern $(BIN_DIR)/size_of_idx 

cleanresults: 
	@echo "Remove result files"
	@rm -f $(TIME_FILES) $(RESULT_FILE) $(INFO_FILES)
	@rm -f $(PATTERNS)

cleanall: clean cleanresults
	@echo "Remove all generated files."
	@rm -f $(TMP_DIR)/* 
	@rm -f $(PAT_DIR)/*
