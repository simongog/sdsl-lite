CC = g++ 
CFLAGS = -O9 -msse4.2 -funroll-loops -fomit-frame-pointer -ffast-math -DNDEBUG
CFLAGS_NOSSE = -O9 -funroll-loops -fomit-frame-pointer -ffast-math -DNDEBUG
CFLAGS_NOOPT = -O9 -DPOPCOUNT_TL -funroll-loops -fomit-frame-pointer -ffast-math -DNDEBUG
LIBS = -lsdsl -ldivsufsort -ldivsufsort64
LIB_DIR = @CMAKE_INSTALL_PREFIX@/lib
INC_DIR = @CMAKE_INSTALL_PREFIX@/include
SRC_DIR = src
TEST_CASES = $(shell find data -name "*[0-9][0-9]MB")

# List of index names. See file csa_typedefs for index type definitions.
INDEXES = FM_HUFF FM_HUFF_RRR15 FM_HUFF_RRR63 FM_HUFF_RRR127 FM_HUFF_RRR255\
		  CSA_SADA FM_RLMN
# List of sample rates for the locate experiment.		  
SAMPLES = 4 8 16 32 128 256
TMP_DIR = tmp

COUNT_RESULT_FILE=results/count.txt

COUNT_EXECS = $(foreach INDEX,$(INDEXES),bin/count_queries_$(INDEX))
COUNT_EXECS_NOSSE = $(foreach INDEX,$(INDEXES),bin/NOSSE/count_queries_$(INDEX))
COUNT_EXECS_HP = $(foreach INDEX,$(INDEXES),bin/HP/count_queries_$(INDEX))
COUNT_EXECS_NOOPT = $(foreach INDEX,$(INDEXES),bin/NOOPT/count_queries_$(INDEX))
BUILD_COUNT_EXECS = $(foreach INDEX,$(INDEXES),bin/build_count_$(INDEX))
GET_INDEX_STRUCTURE_EXECS = $(foreach INDEX,$(INDEXES),bin/get_index_structure_$(INDEX))
TEST_CASES_PATTERNS = $(foreach TEST_CASE,$(TEST_CASES),$(subst data/,pattern/,$(TEST_CASE)).pattern)
#LOCATE_EXECS = $(foreach INDEX,$(INDEXES),$(foreach SAMPLE,$(SAMPLES),bin/locate_queries_$(INDEX).$(SAMPLE)))
COUNT_INDEXES = $(foreach INDEX,$(INDEXES),$(foreach TEST_CASE,$(TEST_CASES),$(TEST_CASE).$(INDEX)))
STRUCTURE_FILES = $(foreach INDEX,$(INDEXES),$(foreach TEST_CASE,$(TEST_CASES),structure/$(subst data/,,$(TEST_CASE)).$(INDEX).json))
TIME_FILES = $(foreach INDEX,$(INDEXES),$(foreach TEST_CASE,$(TEST_CASES),results/$(subst data/,,$(TEST_CASE)).$(INDEX).count))

all: $(BUILD_COUNT_EXECS) $(COUNT_EXECS) $(COUNT_EXECS_NOSSE) $(COUNT_EXECS_NOOPT) $(COUNT_EXECS_HP) \
	 $(GET_INDEX_STRUCTURE_EXECS) bin/genpatterns $(TEST_CASES_PATTERNS)
# $(LOCATE_EXECS)

structure: $(STRUCTURE_FILES)

indexes: $(COUNT_INDEXES)

pattern: $(TEST_CASES_PATTERNS)

timing: $(COUNT_INDEXES) $(TEST_CASES_PATTERNS) $(TIME_FILES)
	cat $(TIME_FILES) > $(COUNT_RESULT_FILE)

results/%.count: $(BUILD_COUNT_EXECS) $(COUNT_EXECS) $(COUNT_EXECS_NOSSE) \
	             $(COUNT_EXECS_NOOPT) $(COUNT_EXECS_HP) $(TEST_CASES_PATTERNS)
	$(eval SUFFIX:=$(suffix $*))				   
	$(eval LOCAL_TEST_CASE:=$(subst $(SUFFIX),,$*))    
	$(eval LOCAL_INDEX:=$(subst .,,$(SUFFIX)))
	$(shell echo "" > $@)
	echo "# test_case = $(LOCAL_TEST_CASE)" >>  $@
	bin/NOOPT/count_queries_$(LOCAL_INDEX) data/$(LOCAL_TEST_CASE) C < pattern/$(LOCAL_TEST_CASE).pattern 2>> $@ 
	echo "# test_case = $(LOCAL_TEST_CASE)" >>  $@
	bin/count_queries_$(LOCAL_INDEX) data/$(LOCAL_TEST_CASE) C < pattern/$(LOCAL_TEST_CASE).pattern 2>> $@ 
	echo "# test_case = $(LOCAL_TEST_CASE)" >>  $@
	bin/NOSSE/count_queries_$(LOCAL_INDEX) data/$(LOCAL_TEST_CASE) C < pattern/$(LOCAL_TEST_CASE).pattern 2>> $@ 
	echo "# test_case = $(LOCAL_TEST_CASE)" >>  $@
	bin/HP/count_queries_$(LOCAL_INDEX) data/$(LOCAL_TEST_CASE) C < pattern/$(LOCAL_TEST_CASE).pattern 2>> $@ 
 

data/%: $(BUILD_COUNT_EXECS)
	$(eval SUFFIX:=$(suffix $*))				   
	$(eval LOCAL_TEST_CASE:=$(subst $(SUFFIX),,$*))    
	$(eval LOCAL_INDEX:=$(subst .,,$(SUFFIX)))  
	bin/build_count_$(LOCAL_INDEX) data/$(LOCAL_TEST_CASE) $(TMP_DIR)/

# Targets for the locate experiment.
structure/%.json: indexes
	@echo $*
	$(eval SUFFIX:=$(suffix $*))				   
	$(eval LOCAL_TEST_CASE:=$(subst data/,,$(subst $(SUFFIX),,$*)))    
	@echo $(LOCAL_TEST_CASE)
	$(eval LOCAL_INDEX:=$(subst .,,$(SUFFIX)))  
	@echo $(LOCAL_INDEX)
	bin/get_index_structure_$(LOCAL_INDEX) data/$(LOCAL_TEST_CASE) > $@ 

pattern/%.pattern: bin/genpatterns
	bin/genpatterns data/$* 20 50000 $@

bin/genpatterns: ${SRC_DIR}/genpatterns.c
	gcc -O3 -o $@ ${SRC_DIR}/genpatterns.c

# Targets for the executables which build the indexes for the count experiments.
bin/build_count_%: $(SRC_DIR)/build_index_sdsl.cpp $(SRC_DIR)/csa_typedefs.hpp
	$(CC) $(CFLAGS) \
					-D$* \
					-L$(LIB_DIR) $(SRC_DIR)/build_index_sdsl.cpp \
					-I$(INC_DIR) \
					-o $@ \
					$(LIBS)

# Targets for the count experiment.
bin/count_queries_%: $(SRC_DIR)/run_queries_sdsl.cpp $(SRC_DIR)/csa_typedefs.hpp 
	$(CC) $(CFLAGS) \
					-D$* \
					-L$(LIB_DIR) $(SRC_DIR)/run_queries_sdsl.cpp \
					-I$(INC_DIR) \
					-o $@ \
					$(LIBS)

# Targets for the count experiment without enabled SSE.
bin/NOSSE/count_queries_%: $(SRC_DIR)/run_queries_sdsl.cpp $(SRC_DIR)/csa_typedefs.hpp 
	$(CC) $(CFLAGS_NOSSE) \
					-D$* \
					-L$(LIB_DIR) $(SRC_DIR)/run_queries_sdsl.cpp \
					-I$(INC_DIR) \
					-o $@ \
					$(LIBS)

# Targets for the count experiment without enabled SSE.
bin/NOOPT/count_queries_%: $(SRC_DIR)/run_queries_sdsl.cpp $(SRC_DIR)/csa_typedefs.hpp 
	$(CC) $(CFLAGS_NOOPT) \
					-D$* \
					-L$(LIB_DIR) $(SRC_DIR)/run_queries_sdsl.cpp \
					-I$(INC_DIR) \
					-o $@ \
					$(LIBS)

# Targets for the count experiment without enabled SSE.
bin/HP/count_queries_%: $(SRC_DIR)/run_queries_sdsl.cpp $(SRC_DIR)/csa_typedefs.hpp 
	$(CC) $(CFLAGS) -DUSE_HP \
					-D$* \
					-L$(LIB_DIR) $(SRC_DIR)/run_queries_sdsl.cpp \
					-I$(INC_DIR) \
					-o $@ \
					$(LIBS)

# Targets for the executables which output the indexes structure.
bin/get_index_structure_%: $(SRC_DIR)/get_index_structure.cpp $(SRC_DIR)/csa_typedefs.hpp 
	$(CC) $(CFLAGS) \
					-D$* \
					-L$(LIB_DIR) $(SRC_DIR)/get_index_structure.cpp \
					-I$(INC_DIR) \
					-o $@ \
					$(LIBS)

# Targets for the locate experiment.
bin/locate_queries_%: $(SRC_DIR)/run_queries_sdsl.cpp $(SRC_DIR)/csa_typedefs.hpp 
	$(eval SUFFIX:=$(suffix $*))				   
	$(eval LOCAL_INDEX:=$(subst $(SUFFIX),,$*))    
	$(eval SA_SAMPLE_RATE:=$(subst .,,$(SUFFIX)))  
	$(CC) $(CFLAGS) \
					-D$(LOCAL_INDEX) \
					-DS_SA=$(SA_SAMPLE_RATE) \
					-DS_ISA=$(SA_SAMPLE_RATE) \
					-L$(LIB_DIR) $(SRC_DIR)/run_queries_sdsl.cpp \
					-I$(INC_DIR) \
					-o $@ \
					$(LIBS)

clean:
	rm -f $(COUNT_EXECS) $(COUNT_EXECS_NOSSE) $(COUNT_EXECS_HP) $(LOCATE_EXECS) \
	   	  $(BUILD_COUNT_EXECS) $(GET_INDEX_STRUCTURE_EXECS) $(COUNT_EXECS_NOOPT) \
		  bin/genpatterns

cleanresults: 
	rm -f $(TIME_FILES) $(COUNT_RESULT_FILE)

cleanall: clean cleanresults
	rm -f $(COUNT_INDEXES) $(STRUCTURE_FILES) $(TEST_CASES_PATTERNS)
	rm -f tmp/*
