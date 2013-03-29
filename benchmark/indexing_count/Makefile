CC = g++ 
CFLAGS = -O9 -funroll-loops -fomit-frame-pointer -ffast-math -DNDEBUG
LIBS = -lsdsl -ldivsufsort -ldivsufsort64
LIB_DIR = /Users/sgog/lib
INC_DIR = /Users/sgog/include
SRC_DIR = src
TMP_DIR = ../tmp
PAT_DIR = pattern
BIN_DIR = bin

# Returns $1-th .-separated part of string $2.
dim = $(word $1, $(subst ., ,$2))

# Returns value stored in column $3 for item with ID $2 in 
# config file $1
config_select=$(shell cat $1 | grep "$2;" | cut -f $3 -d';' )

# Get all IDs from a config file $1
config_ids=$(shell cat $1 | grep -v "^\#" | cut -f 1 -d';')

TC_IDS:=$(call config_ids, test_case.config)
IDX_IDS:=$(call config_ids, index.config)
COMPILE_IDS:=$(call config_ids, compile_options.config)

RESULT_FILE=results/all.txt

QUERY_EXECS = $(foreach IDX_ID,$(IDX_IDS),$(foreach COMPILE_ID,$(COMPILE_IDS),$(BIN_DIR)/query_idx_$(IDX_ID).$(COMPILE_ID)))
BUILD_EXECS = $(foreach IDX_ID,$(IDX_IDS),$(BIN_DIR)/build_idx_$(IDX_ID))
INFO_EXECS = $(foreach IDX_ID,$(IDX_IDS),$(BIN_DIR)/info_$(IDX_ID))
PATTERNS = $(foreach TC_ID,$(TC_IDS),$(PAT_DIR)/$(TC_ID).pattern)
INDEXES = $(foreach IDX_ID,$(IDX_IDS),$(foreach TC_ID,indexes/$(TC_IDS),$(TC_ID).$(IDX_ID)))
INFO_FILES = $(foreach IDX_ID,$(IDX_IDS),$(foreach TC_ID,$(TC_IDS),info/$(TC_ID).$(IDX_ID).json))
TIME_FILES = $(foreach IDX_ID,$(IDX_IDS),$(foreach TC_ID,$(TC_IDS),$(foreach COMPILE_ID,$(COMPILE_IDS),results/$(TC_ID).$(IDX_ID).$(COMPILE_ID))))

all: $(BUILD_EXECS) $(QUERY_EXECS) info pattern

info: $(INFO_FILES)

indexes: $(INDEXES)

pattern: $(PATTERNS) $(BIN_DIR)/genpatterns

timing: $(INDEXES) pattern $(TIME_FILES)
	cat $(TIME_FILES) > $(RESULT_FILE)

results/%: $(BUILD_EXECS) $(QUERY_EXECS) $(PATTERNS)
	$(eval TC_ID:=$(call dim, 1, $*)) 
	$(eval IDX_ID:=$(call dim, 2, $*)) 
	$(eval COMPILE_ID:=$(call dim, 3, $*)) 
	$(eval TC:=$(call config_select,test_case.config,$(TC_ID),2))
	echo "# test_case = $(TC_ID)" >>  $@
	$(BIN_DIR)/query_idx_$(IDX_ID).$(COMPILE_ID) indexes/$(TC_ID) C < $(PAT_DIR)/$(TC_ID).pattern 2>> $@ 
 

# indexes/[TC_ID].[IDX_ID]
indexes/%: $(BUILD_EXECS)
	$(eval TC_ID:=$(call dim,1,$*)) 
	$(eval IDX_ID:=$(call dim,2,$*)) 
	$(eval TC:=$(call config_select,test_case.config,$(TC_ID),2))
	$(BIN_DIR)/build_idx_$(IDX_ID) $(TC) $(TMP_DIR) $@

# info/[TC_ID].[IDX_ID]
info/%.json: indexes
	@echo $*
	$(eval SUFFIX:=$(suffix $*))				   
	$(eval LOCAL_TEST_CASE:=$(subst data/,,$(subst $(SUFFIX),,$*)))    
	@echo $(LOCAL_TEST_CASE)
	$(eval LOCAL_IDX:=$(subst .,,$(SUFFIX)))  
	@echo $(LOCAL_IDX)
	$(BIN_DIR)/info_$(LOCAL_IDX) indexes/$(LOCAL_TEST_CASE) > $@ 

$(PAT_DIR)/%.pattern: $(BIN_DIR)/genpatterns
	@echo $*
	$(eval TC:=$(call config_select,test_case.config,$*,2))
	$(BIN_DIR)/genpatterns $(TC) 20 50000 $@

$(BIN_DIR)/genpatterns: $(SRC_DIR)/genpatterns.c
	gcc -O3 -o $@ $(SRC_DIR)/genpatterns.c

# $(BIN_DIR)/build_idx_[IDX_ID]
$(BIN_DIR)/build_idx_%: $(SRC_DIR)/build_index_sdsl.cpp index.config
	$(eval IDX_TYPE:=$(call config_select,index.config,$*,2))
	$(CC) $(CFLAGS) \
					-DSUF=\"$*\" \
					-DCSA_TYPE="$(IDX_TYPE)" \
					-L$(LIB_DIR) $(SRC_DIR)/build_index_sdsl.cpp \
					-I$(INC_DIR) \
					-o $@ \
					$(LIBS)

# Targets for the count experiment. Pattern $(BIN_DIR)/count_queries_[IDX_ID].[COMPILE_ID]
$(BIN_DIR)/query_idx_%: $(SRC_DIR)/run_queries_sdsl.cpp index.config 
	$(eval IDX_ID:=$(call dim,1,$*)) 
	$(eval COMPILE_ID:=$(call dim,2,$*)) 
	$(eval IDX_TYPE:=$(call config_select,index.config,$(IDX_ID),2))
	$(eval COMPILE_OPTIONS:=$(call config_select,compile_options.config,$(COMPILE_ID),2))
	$(CC) $(CFLAGS) $(COMPILE_OPTIONS) \
					-DSUF=\"$*\" \
					-DCSA_TYPE="$(IDX_TYPE)" \
					-L$(LIB_DIR) $(SRC_DIR)/run_queries_sdsl.cpp \
					-I$(INC_DIR) \
					-o $@ \
					$(LIBS)

# Targets for the executables which output the indexes structure.
$(BIN_DIR)/info_%: $(SRC_DIR)/info.cpp index.config 
	$(eval IDX_TYPE:=$(call config_select,index.config,$*,2))
	$(CC) $(CFLAGS) \
					-DSUF=\"$*\" \
					-DCSA_TYPE="$(IDX_TYPE)" \
					-L$(LIB_DIR) $(SRC_DIR)/info.cpp \
					-I$(INC_DIR) \
					-o $@ \
					$(LIBS)

clean:
	rm -f $(QUERY_EXECS) $(LOCATE_EXECS) $(BUILD_EXECS) $(INFO_EXECS) \
		  $(BIN_DIR)/genpatterns

cleanresults: 
	rm -f $(TIME_FILES) $(RESULT_FILE)

cleanall: clean cleanresults
	rm -f $(INDEXES) $(INFO_FILES) $(PATTERNS)
	rm -f $(TMP_DIR)/* 
	rm -f $(PAT_DIR)/*
