include ../Make.helper
LIBS=-lsdsl -ldivsufsort -ldivsufsort64
SRC_DIR=src
TMP_DIR=../tmp

TC_IDS:=$(call config_ids,test_case.config)
TC_IDS_INT:=$(call config_ids,test_case_int.config)
TC_PATHS:=$(call config_column,test_case.config,2)
TC_PATHS_INT:=$(call config_column,test_case_int.config,2)
IDX_IDS:=$(call config_ids,index.config)
IDX_IDS_INT:=$(call config_ids,index_int.config)
COMPILE_IDS:=$(call config_ids,compile_options.config)
DEFAULT_COMPILE_ID:=$(word 1, $(COMPILE_IDS))

BUILD_EXECS           = $(foreach IDX_ID,$(IDX_IDS),bin/build_cst_$(IDX_ID))
BUILD_EXECS_INT       = $(foreach IDX_ID,$(IDX_IDS_INT),bin/build_int_cst_$(IDX_ID))
CST_EXECS            = $(foreach IDX_ID,$(IDX_IDS),\
						  $(foreach COMPILE_ID,$(COMPILE_IDS),bin/time_and_space_cst_$(IDX_ID).$(COMPILE_ID)))
CST_EXECS_INT        = $(foreach IDX_ID,$(IDX_IDS_INT),\
						  $(foreach COMPILE_ID,$(COMPILE_IDS),bin/time_and_space_int_cst_$(IDX_ID).$(COMPILE_ID)))
NODE_DEPTH_EXECS      = bin/depth_stats bin/depth_stats_int
INDEXES               = $(foreach TC_ID,$(TC_IDS),\
						  $(foreach IDX_ID,$(IDX_IDS),indexes/$(TC_ID).$(IDX_ID).byte))
INDEXES_INT           = $(foreach TC_ID,$(TC_IDS_INT),\
						  $(foreach IDX_ID,$(IDX_IDS_INT),indexes/$(TC_ID).$(IDX_ID).int))
INDEX_FILES           = $(foreach IDX,$(INDEXES),$(IDX))
INDEX_FILES_INT       = $(foreach IDX,$(INDEXES_INT),$(IDX))
RES_FILES             = $(foreach TC_ID,$(TC_IDS),\
						  $(foreach IDX_ID,$(IDX_IDS),\
							$(foreach COMPILE_ID,$(COMPILE_IDS),results/$(TC_ID).$(IDX_ID).$(COMPILE_ID).byte)))
RES_FILES_INT         = $(foreach TC_ID,$(TC_IDS_INT),\
						  $(foreach IDX_ID,$(IDX_IDS_INT),\
							$(foreach COMPILE_ID,$(COMPILE_IDS),results/$(TC_ID).$(IDX_ID).$(COMPILE_ID).int)))
NODE_STATS            = $(foreach TC_ID,$(TC_IDS),stats/$(TC_ID).byte.depth.txt)
NODE_STATS_INT        = $(foreach TC_ID,$(TC_IDS_INT),stats/$(TC_ID).int.depth.txt)

RES_FILE = results/all.txt

all: execs

timing: input\
		$(BUILD_EXECS) $(BUILD_EXECS_INT)\
		$(CST_EXECS) $(CST_EXECS_INT)\
		indexes\
		$(RES_FILES) $(RES_FILES_INT)
	cat $(RES_FILES) $(RES_FILES_INT) > $(RES_FILE)
	@cd visualize; make

depth: input $(NODE_DEPTH_EXECS) $(NODE_STATS) $(NODE_STATS_INT)

input: $(TC_PATHS) $(TC_PATHS_INT)

execs: $(BUILD_EXECS) $(BUILD_EXECS_INT) $(CST_EXECS) $(CST_EXECS_INT) $(NODE_DEPTH_EXECS)

indexes: $(INDEXES) $(INDEXES_INT)

# bin/build_cst_[IDX_ID]
bin/build_cst_%: $(SRC_DIR)/build_cst.cpp
	$(eval IDX_ID:=$(call dim,1,$*))
	$(eval IDX_TYPE:=$(call config_select,index.config,$(IDX_ID),2))
	$(eval COMPILE_ID:=$(DEFAULT_COMPILE_ID))
	$(eval COMPILE_OPTIONS:=$(call config_select,compile_options.config,$(COMPILE_ID),2))
	$(MY_CXX) $(MY_CXX_FLAGS) $(COMPILE_OPTIONS) -DCST_TYPE="$(IDX_TYPE)" -L$(LIB_DIR) -I$(INC_DIR) -o $@ $< $(LIBS)

# bin/build_int_cst_[IDX_ID]
bin/build_int_cst_%: $(SRC_DIR)/build_cst.cpp
	$(eval IDX_ID:=$(call dim,1,$*))
	$(eval IDX_TYPE:=$(call config_select,index_int.config,$(IDX_ID),2))
	$(eval COMPILE_ID:=$(DEFAULT_COMPILE_ID))
	$(eval COMPILE_OPTIONS:=$(call config_select,compile_options.config,$(COMPILE_ID),2))
	$(MY_CXX) $(MY_CXX_FLAGS) $(COMPILE_OPTIONS) -DCST_TYPE="$(IDX_TYPE)" -L$(LIB_DIR) -I$(INC_DIR) -o $@ $< $(LIBS)

# bin/time_and_space_cst_[INDEX_ID].[COMPILE_ID]
bin/time_and_space_cst_%: $(SRC_DIR)/time_and_space_cst.cpp
	$(eval IDX_ID:=$(call dim,1,$*))
	$(eval IDX_TYPE:=$(call config_select,index.config,$(IDX_ID),2))
	$(eval IDX_REF_TYPE:=$(REF_CST_TYPE))
	$(eval COMPILE_ID:=$(call dim,2,$*))
	$(eval COMPILE_OPTIONS:=$(call config_select,compile_options.config,$(COMPILE_ID),2))
	$(MY_CXX) $(MY_CXX_FLAGS) $(COMPILE_OPTIONS) -DCST_TYPE="$(IDX_TYPE)" -L$(LIB_DIR) -I$(INC_DIR) -o $@ $< $(LIBS)

# bin/time_and_space_int_cst_[INDEX_ID].[COMPILE_ID]
bin/time_and_space_int_cst_%: $(SRC_DIR)/time_and_space_cst.cpp
	$(eval IDX_ID:=$(call dim,1,$*))
	$(eval IDX_TYPE:=$(call config_select,index_int.config,$(IDX_ID),2))
	$(eval IDX_REF_TYPE:=$(REF_CST_TYPE_INT))
	$(eval COMPILE_ID:=$(call dim,2,$*))
	$(eval COMPILE_OPTIONS:=$(call config_select,compile_options.config,$(COMPILE_ID),2))
	$(MY_CXX) $(MY_CXX_FLAGS) $(COMPILE_OPTIONS) -DCST_TYPE="$(IDX_TYPE)" -L$(LIB_DIR) -I$(INC_DIR) -o $@ $< $(LIBS)

bin/depth_stats: $(SRC_DIR)/depth_stats.cpp
	$(eval IDX_TYPE:=$(REF_CST_TYPE))
	$(eval COMPILE_ID:=$(DEFAULT_COMPILE_ID))
	$(eval COMPILE_OPTIONS:=$(call config_select,compile_options.config,$(COMPILE_ID),2))
	$(MY_CXX) $(MY_CXX_FLAGS) $(COMPILE_OPTIONS) -DCST_TYPE="$(IDX_TYPE)" -L$(LIB_DIR) -I$(INC_DIR) -o $@ $< $(LIBS)

bin/depth_stats_int: $(SRC_DIR)/depth_stats.cpp
	$(eval IDX_TYPE:=$(REF_CST_TYPE_INT))
	$(eval COMPILE_ID:=$(DEFAULT_COMPILE_ID))
	$(eval COMPILE_OPTIONS:=$(call config_select,compile_options.config,$(COMPILE_ID),2))
	$(MY_CXX) $(MY_CXX_FLAGS) $(COMPILE_OPTIONS) -DCST_TYPE="$(IDX_TYPE)" -L$(LIB_DIR) -I$(INC_DIR) -o $@ $< $(LIBS)

# indexes/[TC_ID].[IDX_ID].byte
indexes/%.byte:
	$(eval TC_ID:=$(call dim,1,$*))
	$(eval IDX_ID:=$(call dim,2,$*))
	$(eval TC_PATH:=$(call config_select,test_case.config,$(TC_ID),2))
	$(eval TIME_FILE:=stats/$(TC_ID).$(IDX_ID).construction.html)
	$(eval STATS_FILE:=stats/$(TC_ID).$(IDX_ID).structure.html)
	@echo "Building index $(IDX_ID) on $(TC_ID) (byte alphabet)"
	@bin/build_cst_$(IDX_ID) $(TC_PATH) $(TMP_DIR) $@ $(TIME_FILE) $(STATS_FILE)

# indexes/[TC_ID].[IDX_ID].int
indexes/%.int:
	$(eval TC_ID:=$(call dim,1,$*))
	$(eval IDX_ID:=$(call dim,2,$*))
	$(eval TC_PATH:=$(call config_select,test_case_int.config,$(TC_ID),2))
	$(eval TIME_FILE:=stats/$(TC_ID).$(IDX_ID).construction.html)
	$(eval STATS_FILE:=stats/$(TC_ID).$(IDX_ID).structure.html)
	@echo "Building index $(IDX_ID) on $(TC_ID) (int alphabet)"
	@bin/build_int_cst_$(IDX_ID) $(TC_PATH) $(TMP_DIR) $@ $(TIME_FILE) $(STATS_FILE)

# Format: results/[TC_ID].[INDEX_ID].[COMPILE_ID].byte
results/%.byte:
	$(eval TC_ID:=$(call dim,1,$*))
	$(eval IDX_ID:=$(call dim,2,$*))
	$(eval COMPILE_ID:=$(call dim,3,$*))
	@echo "Running bin/time_and_space_cst_$(IDX_ID).$(COMPILE_ID) on $(TC_ID)"
	@echo "# TC_ID = $(TC_ID)" >> $@
	@echo "# IDX_ID = $(IDX_ID)" >> $@
	@echo "# COMPILE_ID = $(COMPILE_ID)" >> $@
	@bin/time_and_space_cst_$(IDX_ID).$(COMPILE_ID) indexes/$(TC_ID).$(IDX_ID).byte >> $@

# Format: results/[TC_ID].[INDEX_ID].[COMPILE_ID].int
results/%.int:
	$(eval TC_ID:=$(call dim,1,$*))
	$(eval IDX_ID:=$(call dim,2,$*))
	$(eval COMPILE_ID:=$(call dim,3,$*))
	@echo "Running bin/time_and_space_cst_$(IDX_ID).$(COMPILE_ID) on $(TC_ID)"
	@echo "# TC_ID = $(TC_ID)" >> $@
	@echo "# IDX_ID = $(IDX_ID)" >> $@
	@echo "# COMPILE_ID = $(COMPILE_ID)" >> $@
	@bin/time_and_space_int_cst_$(IDX_ID).$(COMPILE_ID) indexes/$(TC_ID).$(IDX_ID).int >> $@

# stats/[TC_ID].byte.depth.txt
stats/%.byte.depth.txt:
	$(eval TC_ID:=$*)
	$(eval TC_PATH:=$(call config_select,test_case.config,$(TC_ID),2))
	@echo "Generating node depth statistics for $(TC_ID) (byte alphabet)"
	@bin/depth_stats $(TC_PATH) $(TMP_DIR) >> $@

# stats/[TC_ID].int.depth.txt
stats/%.int.depth.txt:
	$(eval TC_ID:=$*)
	$(eval TC_PATH:=$(call config_select,test_case_int.config,$(TC_ID),2))
	@echo "Generating node depth statistics for $(TC_ID) (int alphabet)"
	@bin/depth_stats_int $(TC_PATH) $(TMP_DIR) >> $@

../data/%:
	$(eval URL:=$(call config_filter,test_case.config,$@,4))
	@$(if $(URL),,\
		$(eval URL:=$(call config_filter,test_case_int.config,$@,4))\
		@$(if $(URL),,\
			$(error "No download link nor generation program specified for test case $@") ) )
	@echo "Download input from $(URL) using curl"
	$(eval DEST_DIR:=$(shell dirname $@))
	cd $(DEST_DIR); curl -O $(URL)
	$(eval FILE:=$(DEST_DIR)/$(notdir $(URL)))
	@$(if $(filter-out ".gz",$(FILE)),\
		echo "Extract file $(FILE) using gunzip";\
		gunzip $(FILE))

clean:
	rm -f $(BUILD_EXECS) $(BUILD_EXECS_INT) $(FCST_EXECS) $(FCST_EXECS_INT) $(FCST_VERIFY_EXECS) $(FCST_VERIFY_EXECS_INT) $(TMP_DIR)/*

clean_indexes:
	rm -f $(INDEX_FILES) $(INDEX_FILES_INT)

clean_results:
	rm -f $(RES_FILES) $(RES_FILES_INT)

cleanall: clean clean_indexes clean_results
