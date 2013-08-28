include ../Make.helper
CXX_FLAGS = $(MY_CXX_FLAGS) # in compile_options.config
LIBS = -lsdsl 
SRC_DIR = src
TMP_DIR = ../tmp
# Sampling rate for rrr_vector
t=32
# Block sizes
TC_PATHS:=$(call config_column,test_case.config,2)
TC_IDS:=$(call config_ids,test_case.config)
COMPILE_IDS:=$(call config_ids,compile_options.config)
KS:=$(call config_ids,block_size.config)

all: execs 

input: bin/generate_rnd_bitvector $(TC_PATHS)

RRR_EXECS = $(foreach K,$(KS),\
			  $(foreach COMPILE_ID,$(COMPILE_IDS),bin/rrr_time_and_space_$(K).$(COMPILE_ID)))

RES_FILES = $(foreach TC_ID,$(TC_IDS),\
              $(foreach K,$(KS),\
				$(foreach COMPILE_ID,$(COMPILE_IDS),\
					results/$(TC_ID).$(K).$(COMPILE_ID))))
						
RES_FILE=results/all.txt

# Target for the generator programs for the random bitvectors of different densities
bin/generate_rnd_bitvector: ${SRC_DIR}/generate_rnd_bitvector.cpp
	$(MY_CXX) -O3 $(CXX_FLAGS) $(SRC_DIR)/generate_rnd_bitvector.cpp -L$(LIB_DIR) -I$(INC_DIR) -o $@ $(LIBS)

# Targets for the access/rank/select rrr_vector experiment
# Format: bin/rrr_time_and_space_[R].[COMPILE_ID]
bin/rrr_time_and_space_%: $(SRC_DIR)/rrr_time_and_space.cpp 
	$(eval K:=$(call dim,1,$*)) 
	$(eval COMPILE_ID:=$(call dim,2,$*)) 
	$(eval COMPILE_OPTIONS:=$(call config_select,compile_options.config,$(COMPILE_ID),2))
	$(MY_CXX) $(CXX_FLAGS) $(COMPILE_OPTIONS) -DBLOCK_SIZE=$(K) -L$(LIB_DIR) \
		  $(SRC_DIR)/rrr_time_and_space.cpp -I$(INC_DIR) -o $@ $(LIBS)

execs: $(RRR_EXECS) 


timing: input execs $(RES_FILES) 
	cat $(RES_FILES) > $(RES_FILE)
	@cd visualize; make	

# Format: results/[TC_ID].[K].[COMPILE_ID]
results/%: 
	$(eval TC_ID:=$(call dim,1,$*)) 
	$(eval K:=$(call dim,2,$*)) 
	$(eval COMPILE_ID:=$(call dim,3,$*)) 
	$(eval TC_PATH:=$(call config_select,test_case.config,$(TC_ID),2))
	@echo "Running bin/rrr_time_and_space_$(K).$(COMPILE_ID) on $(TC_ID)"
	@echo "# TC_ID = $(TC_ID)" >> $@
	@echo "# K = $(K)" >> $@
	@echo "# COMPILE_ID = $(COMPILE_ID)" >> $@
	@bin/rrr_time_and_space_$(K).$(COMPILE_ID) $(TC_PATH) $(t) >> $@

# Target for generating a random bitvector of size SIZE and density DENSITY
# Format: data/rnd_[DENSITY].[SIZE]
../data/rnd_%: bin/generate_rnd_bitvector
	$(eval DENSITY:=$(call dim,1,$*)) 
	$(eval SIZE:=$(call dim,2,$*)) 
	@echo "Generating bitvector of density $(DENSITY) and size $(SIZE)"
	@bin/generate_rnd_bitvector $(SIZE) $(DENSITY) $@

include ../Make.download

clean-build:
	@echo "Remove executables"
	rm -f $(RRR_EXECS)

clean:
	rm -f $(RRR_EXECS) bin/generate_rnd_bitvector   

clean_results: 
	rm -f $(RES_FILES) 
		

cleanall: clean clean_results
