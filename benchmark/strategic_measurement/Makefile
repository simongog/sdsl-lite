include ../Make.helper
include ../Make.download

CC = g++ 
LIBS = -lsdsl -ldivsufsort -ldivsufsort64
SRC_DIR = src
TMP_DIR = ../tmp
BIN_DIR = bin
CSA_TYPE="csa_wt<wt_huff<bit_vector,rank_support_v5<>,\
select_support_scan<>,select_support_scan<> >,32,32,\
text_order_sa_sampling<sd_vector<> > >"

TC_PATHS:=$(call config_column,test_case.config,2)
TC_IDS:=$(call config_ids,test_case.config)
PAT_IDS:=$(call config_ids,pattern.config)
MS:=8 16 32 64
SKIPS:=0 1000 2000 3000 4000 5000 6000 7000 8000 9000
MAX_PAT:=1000

OCC_FILES = $(foreach TC_ID,$(TC_IDS),results/$(TC_ID).occ)
OCC_FILE  = results/occ.txt
PAT_FILES = $(foreach PAT_ID,$(PAT_IDS),\
              $(foreach TC_ID,$(TC_IDS),pattern/$(TC_ID).$(PAT_ID)))
CAND_FILES = $(foreach PAT_ID,$(PAT_IDS),\
              $(foreach TC_ID,$(TC_IDS),results/$(TC_ID).$(PAT_ID).cand))
CAND_FILE = results/cand.txt
COUNT_RES = $(foreach TC_ID,$(TC_IDS),\
              $(foreach SKIP,$(SKIPS),\
                $(foreach PAT_ID,$(PAT_IDS),results/$(TC_ID).$(PAT_ID).$(SKIP).count)))
COUNT_FILE= results/count.txt
LOCATE_RES = $(foreach TC_ID,$(TC_IDS),\
              $(foreach SKIP,$(SKIPS),\
                $(foreach PAT_ID,$(PAT_IDS),results/$(TC_ID).$(PAT_ID).$(SKIP).locate)))
LOCATE_FILE= results/locate.txt

INDEXES   = $(foreach TC_ID,$(TC_IDS),indexes/$(TC_ID))
LCP_DIST = $(foreach TC_ID,$(TC_IDS),results/lcp-dist-$(TC_ID))

COMP_FILES  = $(addsuffix .z.info,$(TC_PATHS))

lcp-info:
	@echo $(LCP_DIST)

occ-table: input $(OCC_FILES)
	@cat $(OCC_FILES) > $(OCC_FILE)
	@cd visualize; make occ.pdf

all: input $(OCC_FILES) 
	@cat $(OCC_FILES) > $(OCC_FILE)
	@cd visualize; make

compression: input $(COMP_FILES)

input: $(TC_PATHS)

pattern: $(BIN_DIR)/genpatterns $(BIN_DIR)/stratified_pattern input $(PAT_FILES) bin/extract_candidates $(CAND_FILES)
	@cat $(CAND_FILES) > $(CAND_FILE)
	cd visualize; make cand.pdf

#$(CAND_FILE): bin/extract_candidates $(CAND_FILES)
	

timing: count-timing locate-timing

count-timing: $(INDEXES) $(COUNT_RES) 
	@cat $(COUNT_RES) > $(COUNT_FILE)
	cd visualize; make count.pdf

locate-timing: $(INDEXES) $(LOCATE_RES)
	@cat $(LOCATE_RES) > $(LOCATE_FILE)
	cd visualize; make locate.pdf

lcp-dist: input bin/lcp-distribution $(LCP_DIST)

# Format results/lcp-dist-[TC_ID]
results/lcp-dist-%: bin/lcp-distribution 
	$(eval TC_ID:=$*)
	$(eval TC:=$(call config_select,test_case.config,$(TC_ID),2))
	@echo "Calculating lcp distribution for $(TC_ID) ($(TC))"
	@bin/lcp-distribution $(TC) ../tmp > $@

# Format: indexes/[TC_ID]
indexes/%: $(BIN_DIR)/build_index_sdsl
	$(eval TC_ID:=$*)
	$(eval TC:=$(call config_select,test_case.config,$(TC_ID),2))
	@echo "Building index for $(TC)"
	@$(BIN_DIR)/build_index_sdsl $(TC) $(TMP_DIR) $@

# Format: pattern/[TC_ID].[PAT_ID]
pattern/%-random-10000: input $(BIN_DIR)/genpatterns
	$(eval TC_ID:=$(call dim,1,$*))
	$(eval PAT_ID:=$(call dim,2,$*)-random-10000)
	$(eval PAT_LEN:=$(call config_select,pattern.config,$(PAT_ID),2))
	$(eval PAT_NR:=$(call config_select,pattern.config,$(PAT_ID),5))
	$(eval TC_PATH:=$(call config_select,test_case.config,$(TC_ID),2))
	$(BIN_DIR)/genpatterns $(TC_PATH) $(PAT_LEN) $(PAT_NR) $@

# Format: pattern/[TC_ID].[PAT_ID]
pattern/%: input $(BIN_DIR)/stratified_pattern
	$(eval TC_ID:=$(call dim,1,$*))
	$(eval PAT_ID:=$(call dim,2,$*))
	$(eval TC_PATH:=$(call config_select,test_case.config,$(TC_ID),2))
	$(eval PAT_LEN:=$(call config_select,pattern.config,$(PAT_ID),2))
	$(eval PAT_MIN_OCC:=$(call config_select,pattern.config,$(PAT_ID),3))
	$(eval PAT_MAX_OCC:=$(call config_select,pattern.config,$(PAT_ID),4))
	$(eval PAT_NR:=$(call config_select,pattern.config,$(PAT_ID),5))
	$(BIN_DIR)/stratified_pattern $(TC_PATH) $(TMP_DIR) $(PAT_LEN) \
                                  $(PAT_MIN_OCC) $(PAT_MAX_OCC) $(PAT_NR) $@

$(BIN_DIR)/genpatterns: $(SRC_DIR)/genpatterns.c
	@echo "Building Pizza&Chili pattern generation program"
	@gcc -O3 -o $@ $(SRC_DIR)/genpatterns.c

$(BIN_DIR)/random: $(SRC_DIR)/random.c
	@echo "Building Moffat's random text generation program"
	@gcc -O3 -o $@ $(SRC_DIR)/random.c

$(BIN_DIR)/%: $(SRC_DIR)/%.cpp
	@echo "Building program $@"
	@$(CC) $(CFLAGS) -O3 -DNDEBUG -DCSA_TYPE=$(CSA_TYPE) -msse4.2 -L$(LIB_DIR) -I$(INC_DIR) $< -o $@ $(LIBS)

data/rnd.%: bin/random
	$(eval SIGMA:=$(call dim,1,$*))
	$(eval BASELEN:=$(call dim,2,$*))
	$(eval REPEATS:=$(call dim,3,$*))
	@echo "Generate random string with sigma=$(SIGMA) baselen=$(BASELEN) repeats=$(REPEATS)"
	@bin/random $(SIGMA) $(BASELEN) $(REPEATS) > $@

# Format: results/[TC_ID].[PAT_ID].[SKIP].count
results/%.count: $(BIN_DIR)/run_queries_sdsl
	$(eval TC_ID:=$(call dim,1,$*))
	$(eval PAT_ID:=$(call dim,2,$*))
	$(eval SKIP:=$(call dim,3,$*))
	$(eval PAT_M:=$(call config_select,pattern.config,$(PAT_ID),2))
	$(eval PAT_MIN_OCC:=$(call config_select,pattern.config,$(PAT_ID),3))
	$(eval PAT_MAX_OCC:=$(call config_select,pattern.config,$(PAT_ID),4))
	$(eval PAT_NR:=$(call config_select,pattern.config,$(PAT_ID),5))
	$(eval TC_PATH:=$(call config_select,test_case.config,$(TC_ID),2))
	@echo "Use patterns $(PAT_ID) to run count queries on $(TC_ID)"
	@echo "# TC_ID = $(TC_ID)" >> $@
	@echo "# PAT_M = $(PAT_M)" >> $@
	@echo "# PAT_MIN_OCC = $(PAT_MIN_OCC)" >> $@
	@echo "# PAT_MAX_OCC = $(PAT_MAX_OCC)" >> $@
	@echo "# PAT_NR = $(PAT_NR)" >> $@
	@$(BIN_DIR)/run_queries_sdsl \
	indexes/$(TC_ID) C $(MAX_PAT) $(SKIP) < pattern/$(TC_ID).$(PAT_ID) 2>> $@ 
 
# Format: results/[TC_ID].[PAT_ID].[SKIP].locate
results/%.locate: $(BIN_DIR)/run_queries_sdsl
	$(eval TC_ID:=$(call dim,1,$*))
	$(eval PAT_ID:=$(call dim,2,$*))
	$(eval SKIP:=$(call dim,3,$*))
	$(eval PAT_M:=$(call config_select,pattern.config,$(PAT_ID),2))
	$(eval PAT_MIN_OCC:=$(call config_select,pattern.config,$(PAT_ID),3))
	$(eval PAT_MAX_OCC:=$(call config_select,pattern.config,$(PAT_ID),4))
	$(eval TC_PATH:=$(call config_select,test_case.config,$(TC_ID),2))
	@echo "Use patterns $(PAT_ID) to run locate queries on $(TC_ID)"
	@echo "# TC_ID = $(TC_ID)" >> $@
	@echo "# PAT_M = $(PAT_M)" >> $@
	@echo "# PAT_MIN_OCC = $(PAT_MIN_OCC)" >> $@
	@echo "# PAT_MAX_OCC = $(PAT_MAX_OCC)" >> $@
	@$(BIN_DIR)/run_queries_sdsl \
	indexes/$(TC_ID) L $(MAX_PAT) $(SKIP) < pattern/$(TC_ID).$(PAT_ID) 2>> $@ > $@.info

# Format: results/[TC_ID].[PAT_ID].cand
results/%.cand: $(BIN_DIR)/extract_candidates
	$(eval TC_ID:=$(call dim,1,$*))
	$(eval PAT_ID:=$(call dim,2,$*))
	$(eval PAT_M:=$(call config_select,pattern.config,$(PAT_ID),2))
	$(eval PAT_MIN_OCC:=$(call config_select,pattern.config,$(PAT_ID),3))
	$(eval PAT_MAX_OCC:=$(call config_select,pattern.config,$(PAT_ID),4))
	@echo "Extract candidates for $(PAT_ID) on $(TC_ID)"
	@echo "# TC_ID = $(TC_ID)" >> $@
	@echo "# PAT_M = $(PAT_M)" >> $@
	@echo "# PAT_MIN_OCC = $(PAT_MIN_OCC)" >> $@
	@echo "# PAT_MAX_OCC = $(PAT_MAX_OCC)" >> $@
	$(eval CAND:=$(shell $(BIN_DIR)/extract_candidates < pattern/$(TC_ID).$(PAT_ID) ))
	@echo "# CANDIDATES = $(CAND)" >> $@

# Format: results/[TC_ID].occ
results/%.occ: $(BIN_DIR)/occurrence_statistics
	$(eval TC_ID:=$*) 
	$(eval TC_PATH:=$(call config_select,test_case.config,$(TC_ID),2))
	@echo "Processing $(TC_ID) with pattern length $(MS)"
	@echo $(MS) | ${BIN_DIR}/occurrence_statistics $(TC_PATH) $(TMP_DIR) $(TC_ID) > $@
