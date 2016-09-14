include ../Make.helper
CXX_FLAGS = $(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(INC_DIR) -L$(LIB_DIR)  
LIBS = -lsdsl -ldivsufsort -ldivsufsort64
SRC_DIR = src
TMP_DIR = ../tmp
PAT_DIR = pattern
BIN_DIR = bin

TC_PATHS:=$(call config_column,test_case.config,2)
TC_PATHS_INT:=$(call config_column,test_case_int.config,2)
TC_DICT_PATHS:=$(call config_column,dic.config,2)
TC_IDS:=$(call config_ids,test_case.config)
TC_IDS_INT:=$(call config_ids,test_case_int.config)
IDX_IDS:=$(call config_ids,index.config)
IDX_IDS_INT:=$(call config_ids,index_int.config)
PAT_LENS:=$(call config_column,pattern_length.config,1)
PAT_LENS_INT:=$(call config_column,pattern_length_int.config,1)

RESULT_FILE=results/all.txt
RESULT_FILE_INT=results/all_int.txt

QUERY_EXECS     = $(foreach IDX_ID,$(IDX_IDS),$(BIN_DIR)/query_idx_$(IDX_ID))
QUERY_EXECS_INT = $(foreach IDX_ID,$(IDX_IDS_INT),$(BIN_DIR)/query_int_idx_$(IDX_ID))
BUILD_EXECS     = $(foreach IDX_ID,$(IDX_IDS),$(BIN_DIR)/build_idx_$(IDX_ID))
BUILD_EXECS_INT = $(foreach IDX_ID,$(IDX_IDS_INT),$(BIN_DIR)/build_int_idx_$(IDX_ID))
SIZES_EXECS     = $(foreach IDX_ID,$(IDX_IDS),$(BIN_DIR)/size_of_idx_$(IDX_ID))
SIZES_EXECS_INT = $(foreach IDX_ID,$(IDX_IDS_INT),$(BIN_DIR)/size_of_int_idx_$(IDX_ID))
PATTERNS        = $(foreach TC_ID,$(TC_IDS),\
				    $(foreach PAT_LEN,$(PAT_LENS),$(PAT_DIR)/$(TC_ID).$(PAT_LEN).pattern))
PATTERNS_INT    = $(foreach TC_ID,$(TC_IDS_INT),\
				    $(foreach PAT_LEN,$(PAT_LENS_INT),$(PAT_DIR)/$(TC_ID).$(PAT_LEN).pattern.int))
INDEXES         = $(foreach IDX_ID,$(IDX_IDS),\
		            $(foreach TC_ID,$(TC_IDS),indexes/$(TC_ID).$(IDX_ID).byte))
INDEXES_INT     = $(foreach IDX_ID,$(IDX_IDS_INT),\
		            $(foreach TC_ID,$(TC_IDS_INT),indexes/$(TC_ID).$(IDX_ID).int))
SIZES           = $(foreach IDX_ID,$(IDX_IDS),\
		            $(foreach TC_ID,$(TC_IDS),info/$(TC_ID).$(IDX_ID).size))
SIZES_INT       = $(foreach IDX_ID,$(IDX_IDS_INT),\
		            $(foreach TC_ID,$(TC_IDS_INT),info/$(TC_ID).$(IDX_ID).size.int))
HTML            = $(foreach IDX_ID,$(IDX_IDS),\
		            $(foreach TC_ID,$(TC_IDS),info/$(TC_ID).$(IDX_ID).byte.html))
HTML_INT        = $(foreach IDX_ID,$(IDX_IDS_INT),\
		            $(foreach TC_ID,$(TC_IDS_INT),info/$(TC_ID).$(IDX_ID).int.html))
TIME_FILES      = $(foreach IDX_ID,$(IDX_IDS),\
		            $(foreach TC_ID,$(TC_IDS),\
					  $(foreach PAT_LEN,$(PAT_LENS),results/$(TC_ID).$(IDX_ID).$(PAT_LEN).byte)))
TIME_FILES_INT  = $(foreach IDX_ID,$(IDX_IDS_INT),\
		            $(foreach TC_ID,$(TC_IDS_INT),\
					  $(foreach PAT_LEN,$(PAT_LENS_INT),results/$(TC_ID).$(IDX_ID).$(PAT_LEN).int)))
COMP_FILES      = $(addsuffix .z.info,$(TC_PATHS) $(TC_PATHS_INT))

HELPER_BINS     = $(BIN_DIR)/gen_pattern $(BIN_DIR)/gen_pattern_int \
		          $(BIN_DIR)/word_pat2char_pat

all: $(BUILD_EXECS) $(BUILD_EXECS_INT) \
	 $(QUERY_EXECS) $(QUERY_EXECS_INT) \
	 $(SIZES_EXECS) $(SIZES_EXECS_INT) \
	 $(HELPER_BINS)

info: $(SIZES_EXECS) $(SIZES_EXECS_INT) $(SIZES) \
	  $(SIZES_INT) $(HTML_INT) $(HTML) \
	  info/sizes.txt \
	  info/sizes_int.txt

info/sizes.txt: $(SIZES)
	@cat $(SIZES) > $@

info/sizes_int.txt: $(SIZES)
	@cat $(SIZES_INT) > $@

indexes: $(INDEXES) $(INDEXES_INT)

input: $(TC_PATHS) $(TC_PATHS_INT) $(TC_DICT_PATHS)

pattern: input $(PATTERNS) $(BIN_DIR)/gen_pattern $(PATTERNS_INT) $(BIN_DIR)/gen_pattern_int

compression: input $(COMP_FILES)

timing: input indexes pattern $(TIME_FILES) $(TIME_FILES_INT) compression info
	@cat $(TIME_FILES) > $(RESULT_FILE)
	@cat $(TIME_FILES_INT) > $(RESULT_FILE_INT)
	@cd visualize; make

# results/[TC_ID].[IDX_ID].[PAT_LEN].byte
results/%.byte: $(BUILD_EXECS) $(QUERY_EXECS) $(INDEXES) $(PATTERNS)
	$(eval TC_ID:=$(call dim,1,$*)) 
	$(eval IDX_ID:=$(call dim,2,$*)) 
	$(eval PAT_LEN:=$(call dim,3,$*)) 
	$(eval TC_NAME:=$(call config_select,test_case.config,$(TC_ID),3))
	@echo "# TC_ID = $(TC_ID)" > $@
	@echo "# IDX_ID = $(IDX_ID)" >> $@
	@echo "# test_case = $(TC_NAME)" >>  $@
	@echo "Run timing for $(IDX_ID) on $(TC_ID) with patterns of length $(PAT_LEN)"
	@$(BIN_DIR)/query_idx_$(IDX_ID) indexes/$(TC_ID).$(IDX_ID).byte \
		     $(PAT_DIR)/$(TC_ID).$(PAT_LEN).pattern >> $@ 

# results/[TC_ID].[IDX_ID].[PAT_LEN].byte
results/%.int: $(BUILD_EXECS) $(QUERY_EXECS_INT) $(INDEXES_INT) $(PATTERNS_INT)
	$(eval TC_ID:=$(call dim,1,$*)) 
	$(eval IDX_ID:=$(call dim,2,$*)) 
	$(eval PAT_LEN:=$(call dim,3,$*)) 
	$(eval TC_NAME:=$(call config_select,test_case_int.config,$(TC_ID),3))
	@echo "# TC_ID = $(TC_ID)" > $@
	@echo "# IDX_ID = $(IDX_ID)" >> $@
	@echo "# test_case = $(TC_NAME)" >>  $@
	@echo "Run timing for $(IDX_ID) on $(TC_ID) with patterns of length $(PAT_LEN)"
	@$(BIN_DIR)/query_int_idx_$(IDX_ID) indexes/$(TC_ID).$(IDX_ID).int \
		     $(PAT_DIR)/$(TC_ID).$(PAT_LEN).pattern.int >> $@ 
 
 

# indexes/[TC_ID].[IDX_ID].byte
indexes/%.byte: $(BUILD_EXECS)
	$(eval TC_ID:=$(call dim,1,$*)) 
	$(eval IDX_ID:=$(call dim,2,$*)) 
	$(eval TC:=$(call config_select,test_case.config,$(TC_ID),2))
	@echo "Building index $(IDX_ID) on $(TC)"
	@$(BIN_DIR)/build_idx_$(IDX_ID) $(TC) $(TMP_DIR) $@

# indexes/[TC_ID].[IDX_ID].int
indexes/%.int: $(BUILD_EXECS_INT)
	$(eval TC_ID:=$(call dim,1,$*)) 
	$(eval IDX_ID:=$(call dim,2,$*)) 
	$(eval TC:=$(call config_select,test_case_int.config,$(TC_ID),2))
	@echo "Building index $(IDX_ID) on $(TC)"
	@$(BIN_DIR)/build_int_idx_$(IDX_ID) $(TC) $(TMP_DIR) $@

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
	@$(BIN_DIR)/size_of_idx_$(IDX_ID) indexes/$(TC_ID).$(IDX_ID).byte >> $@

# info/[TC_ID].[IDX_ID]
info/%.size.int: $(INDEXES_INT)
	$(eval TC_ID:=$(call dim,1,$*)) 
	$(eval IDX_ID:=$(call dim,2,$*)) 
	$(eval TC:=$(call config_select,test_case_int.config,$(TC_ID),2))
	$(eval SIZE:=$(call file_size,$(TC)))
	@echo "# TC_ID = $(TC_ID)" > $@
	@echo "# IDX_ID = $(IDX_ID)" >> $@
	@echo "# text_size = $(SIZE)" >> $@
	@echo "Get size of index for $(IDX_ID) on $(TC_ID)"
	@$(BIN_DIR)/size_of_int_idx_$(IDX_ID) indexes/$(TC_ID).$(IDX_ID).int >> $@

# info/[TC_ID].[IDX_ID]
info/%.byte.html: $(INDEXES)
	$(eval TC_ID:=$(call dim,1,$*)) 
	$(eval IDX_ID:=$(call dim,2,$*)) 
	$(eval TC:=$(call config_select,test_case.config,$(TC_ID),2))
	$(eval SIZE:=$(call file_size,$(TC)))
	@echo "# TC_ID = $(TC_ID)" > $@
	@echo "# IDX_ID = $(IDX_ID)" >> $@
	@echo "# text_size = $(SIZE)" >> $@
	@echo "Get html of index for $(IDX_ID) on $(TC_ID)"
	@$(BIN_DIR)/size_of_idx_$(IDX_ID) indexes/$(TC_ID).$(IDX_ID).byte $@

# info/[TC_ID].[IDX_ID]
info/%.int.html: $(INDEXES)
	$(eval TC_ID:=$(call dim,1,$*)) 
	$(eval IDX_ID:=$(call dim,2,$*)) 
	$(eval TC:=$(call config_select,test_case_int.config,$(TC_ID),2))
	$(eval SIZE:=$(call file_size,$(TC)))
	@echo "# TC_ID = $(TC_ID)" > $@
	@echo "# IDX_ID = $(IDX_ID)" >> $@
	@echo "# text_size = $(SIZE)" >> $@
	@echo "Get html of index for $(IDX_ID) on $(TC_ID)"
	@$(BIN_DIR)/size_of_int_idx_$(IDX_ID) indexes/$(TC_ID).$(IDX_ID).int $@



# $(PAT_DIR)/[TC_ID].[PAT_LEN].pattern.int
$(PAT_DIR)/%.pattern.int: $(BIN_DIR)/gen_pattern_int $(BIN_DIR)/word_pat2char_pat
	@echo $*
	$(eval TC_ID:=$(call dim,1,$*)) 
	$(eval PAT_LEN:=$(call dim,2,$*)) 
	$(eval TC:=$(call config_select,test_case_int.config,$(TC_ID),2))
	$(BIN_DIR)/gen_pattern_int $(TC) $(TMP_DIR)/$(TC_ID).pat.csa $(TMP_DIR) $(PAT_LEN) 200 $@
	$(eval DIC_PATH:=$(call config_select,dic.config,$(TC_ID),2))
	$(BIN_DIR)/word_pat2char_pat $@ $(DIC_PATH) > $@.txt

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

$(BIN_DIR)/gen_pattern_int: $(SRC_DIR)/gen_pattern.cpp 
	@echo "Build pattern generation program"
	$(MY_CXX) $(CXX_FLAGS) $(SRC_DIR)/gen_pattern.cpp \
	          -DINT_ALPHABET \
			  -L$(LIB_DIR) -I$(INC_DIR) -o $@ $(LIBS)

$(BIN_DIR)/size_of_idx_%: $(SRC_DIR)/size_of_idx.cpp
	$(eval IDX_TYPE:=$(call config_select,index.config,$*,2))
	@echo "Build size info program"
	$(MY_CXX) $(CXX_FLAGS) $(SRC_DIR)/size_of_idx.cpp \
				-DIDX_TYPE="$(IDX_TYPE)" \
				-L$(LIB_DIR) -I$(INC_DIR) -o $@ $(LIBS)

$(BIN_DIR)/size_of_int_idx_%: $(SRC_DIR)/size_of_idx.cpp
	$(eval IDX_TYPE:=$(call config_select,index_int.config,$*,2))
	@echo "Build size info program"
	$(MY_CXX) $(CXX_FLAGS) $(SRC_DIR)/size_of_idx.cpp \
				-DIDX_TYPE="$(IDX_TYPE)" \
				-L$(LIB_DIR) -I$(INC_DIR) -o $@ $(LIBS)

# $(BIN_DIR)/build_idx_[IDX_ID]
$(BIN_DIR)/build_idx_%: $(SRC_DIR)/build_idx.cpp index.config
	$(eval IDX_TYPE:=$(call config_select,index.config,$*,2))
	@echo "Compiling build_idx_$*"
	$(MY_CXX) $(CXX_FLAGS) \
					-DIDX_TYPE="$(IDX_TYPE)" \
					-L$(LIB_DIR) $(SRC_DIR)/build_idx.cpp \
					-I$(INC_DIR) -o $@ $(LIBS)

# $(BIN_DIR)/build_idx_[IDX_ID]
$(BIN_DIR)/build_int_idx_%: $(SRC_DIR)/build_idx.cpp index_int.config
	$(eval IDX_TYPE:=$(call config_select,index_int.config,$*,2))
	@echo "Compiling build_int_idx_$*"
	$(MY_CXX) $(CXX_FLAGS) \
					-DIDX_TYPE="$(IDX_TYPE)" \
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

# Targets for the count experiment. $(BIN_DIR)/query_idx_[IDX_ID]
$(BIN_DIR)/query_int_idx_%: $(SRC_DIR)/query_idx.cpp index_int.config 
	$(eval IDX_TYPE:=$(call config_select,index_int.config,$*,2))
	@echo "Compiling query_int_idx_$*"
	$(MY_CXX) $(CXX_FLAGS) \
					-DIDX_ID=\"$*\" -DIDX_TYPE="$(IDX_TYPE)" \
					-L$(LIB_DIR) $(SRC_DIR)/query_idx.cpp \
					-I$(INC_DIR) -o $@ $(LIBS)

$(BIN_DIR)/word_pat2char_pat: $(SRC_DIR)/word_pat2char_pat.cpp
	@echo "Compiling word_pat2char_pat.cpp"
	$(MY_CXX) $(CXX_FLAGS) \
			  $(SRC_DIR)/word_pat2char_pat.cpp \
			  -o $@ 



dic/%: 
	$(eval URL:=$(call config_filter,dic.config,$@,4))
	@$(if $(URL),,\
		$(error "No download link for dictionary $@") )
	@echo "Download dictionary from $(URL) using curl"
	$(eval DEST_DIR:=$(shell dirname $@))
	cd $(DEST_DIR); curl -O $(URL)
	$(eval FILE:=$(DEST_DIR)/$(notdir $(URL)))
	@$(if $(filter-out ".gz",$(FILE)),\
		echo "Extract file $(FILE) using gunzip";\
		gunzip $(FILE))

../data/%.int.sdsl: 
	$(eval URL:=$(call config_filter,test_case_int.config,$@,4))
	@$(if $(URL),,\
		$(error "No download link nor generation program specified for test case $@") )
	@echo "Download input from $(URL) using curl"
	$(eval DEST_DIR:=$(shell dirname $@))
	cd $(DEST_DIR); curl -O $(URL)
	$(eval FILE:=$(DEST_DIR)/$(notdir $(URL)))
	@$(if $(filter-out ".gz",$(FILE)),\
		echo "Extract file $(FILE) using gunzip";\
		gunzip $(FILE))


include ../Make.download

clean-build:
	@echo "Remove executables"
	rm -f $(QUERY_EXECS) $(SIZES_EXECS) $(BUILD_EXECS) \
		  $(QUERY_EXECS_INT) $(SIZES_EXECS_INT) $(BUILD_EXECS_INT) \
		  $(HELPER_BINS)

clean: clean-build
	@echo "Remove info files and indexes"
	rm -f $(INDEXES) $(INDEXES_INT)

cleanresults: 
	@echo "Remove result files and pattern"
	rm -f $(TIME_FILES) \
		  $(RESULT_FILE) \
		  $(TIME_FILES_INT) \
		  $(RESULT_FILE_INT) \
	      $(SIZES) $(SIZES_INT) \
		  $(HTML) $(HTML_INT)
	rm -f $(PATTERNS) $(PATTERNS_INT)

cleanall: clean cleanresults
	@echo "Remove all generated files."
	rm -f $(TMP_DIR)/* 
	rm -f $(PAT_DIR)/*
	cd visualize; make clean
