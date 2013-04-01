include ../Make.helper
CC = g++ 
CFLAGS = -msse4.2 -O9 -funroll-loops -fomit-frame-pointer -ffast-math -DNDEBUG
LIBS = -lsdsl -ldivsufsort -ldivsufsort64
SRC_DIR = src
TMP_DIR = ../tmp
PAT_DIR = pattern
BIN_DIR = bin

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
PATTERNS    = $(foreach TC_ID,$(TC_IDS),$(PAT_DIR)/$(TC_ID).pattern)
INDEXES     = $(foreach IDX_ID,$(IDX_IDS),\
		        $(foreach TC_ID,$(TC_IDS),\
				  $(foreach SAMPLE_ID,$(SAMPLE_IDS),indexes/$(TC_ID).$(IDX_ID).$(SAMPLE_ID))))
INFO_FILES  = $(foreach IDX_ID,$(IDX_IDS),\
		      	$(foreach TC_ID,$(TC_IDS),\
				  $(foreach SAMPLE_ID,$(SAMPLE_IDS),info/$(TC_ID).$(IDX_ID).$(SAMPLE_ID).json)))
TIME_FILES  = $(foreach IDX_ID,$(IDX_IDS),\
		        $(foreach TC_ID,$(TC_IDS),\
				  $(foreach SAMPLE_ID,$(SAMPLE_IDS),results/$(TC_ID).$(IDX_ID).$(SAMPLE_ID))))
COMP_FILES  = $(foreach TC_ID,$(TC_IDS),results/compression_$(TC_ID))

all: $(BUILD_EXECS) $(QUERY_EXECS) info pattern

info: $(INFO_EXECS) $(INFO_FILES) 

indexes: $(INDEXES)

pattern: $(PATTERNS) $(BIN_DIR)/pattern_random

compression: $(COMP_FILES)

timing: $(INDEXES) pattern $(TIME_FILES) compression
	@cat $(TIME_FILES) > $(RESULT_FILE)
	@cd visualize; make

# results/baseline_[TC_ID]
results/compression_%: test_case.config
	$(eval TC_ID:=$(call dim,1,$*)) 
	$(eval TC_PATH:=$(call config_select,test_case.config,$(TC_ID),2))
	@echo "Get compression of xz for $(TC_ID)"
	$(eval TC_XZ:=$(shell basename $(TC_PATH)).xz)
	$(shell xz -9 -z -k -c $(TC_PATH) > $(TC_XZ))  
	$(eval XZ_SIZE:=$(call file_size,$(TC_XZ))) 
	$(shell rm $(TC_XZ))
	@echo "Get compression of gizp for $(TC_ID)"
	$(eval TC_GZ:=$(shell basename $(TC_PATH)).gz)
	$(shell gzip -9 -c $(TC_PATH) > $(TC_GZ))  
	$(eval GZ_SIZE:=$(call file_size,$(TC_GZ))) 
	$(shell rm $(TC_GZ))
	$(eval SIZE:=$(call file_size,$(TC_PATH)))
	$(eval XZ_RATIO:=$(shell echo "scale=2;100*$(XZ_SIZE)/$(SIZE)" | bc -q))
	$(eval GZ_RATIO:=$(shell echo "scale=2;100*$(GZ_SIZE)/$(SIZE)" | bc -q))
	@echo "xz;$(XZ_RATIO);xz -9\ngzip;$(GZ_RATIO);gzip -9" > $@
	

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
		indexes/$(TC_ID) L < $(PAT_DIR)/$(TC_ID).pattern 2>> $@ 
 

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

# $(PAT_DIR)/[TC_ID].pattern 
$(PAT_DIR)/%.pattern: $(BIN_DIR)/pattern_random
	@echo "Generating pattern for $*"
	$(eval TC:=$(call config_select,test_case.config,$*,2))
	@$(BIN_DIR)/pattern_random $(TC) $(TMP_DIR) 5 2000000 $@ $@.occ

$(BIN_DIR)/pattern_random: $(SRC_DIR)/pattern_random.cpp
	@echo "Build pattern generation program"
	@$(CC) $(CFLAGS) -L$(LIB_DIR) -I$(INC_DIR) \
					 $(SRC_DIR)/pattern_random.cpp  -o $@ $(LIBS)

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
	@$(CC) $(CFLAGS) -DSUF=\"$*\" -DCSA_TYPE="$(IDX_TYPE)" \
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
	@$(CC) $(CFLAGS) -DSUF="$(IDX_ID).$(SAMPLE_ID)" -DCSA_TYPE="$(IDX_TYPE)" \
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
	@$(CC) $(CFLAGS) -DSUF=\"$*\" -DCSA_TYPE="$(IDX_TYPE)" \
					-L$(LIB_DIR) $(SRC_DIR)/info.cpp \
					-I$(INC_DIR) -o $@ $(LIBS)

clean:
	@echo "Remove executables"
	@rm -f $(QUERY_EXECS) $(LOCATE_EXECS) $(BUILD_EXECS) $(INFO_EXECS) \
		   $(BIN_DIR)/pattern_random

cleanresults: 
	@echo "Remove result files"
	@rm -f $(TIME_FILES) $(RESULT_FILE)

cleanall: clean cleanresults
	@echo "Remove all generated files."
	@rm -f $(INDEXES) $(INFO_FILES) $(PATTERNS)
	@rm -f $(TMP_DIR)/* 
	@rm -f $(PAT_DIR)/*
	@rm -f $(BIN_DIR)/pattern_random
