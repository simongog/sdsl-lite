include ../Make.helper
include ../Make.download

MEASURE_ENERGY = 1 
MEASURE_CMD = $(if $(MEASURE_ENERGY), likwid-perfctr -C S0:0 -g ENERGY -m,)

CXX_FLAGS = $(MY_CXX_FLAGS) # in compile_options.config
LIBS = -lsdsl 
SRC_DIR = src
TMP_DIR = ../tmp

.SECONDARY:

TC_IDS:=$(call config_ids,test_case.config)
TC_COLLECTIONS:=#$(foreach TC_ID,$(TC_IDS),collections/$(TC_ID)/text.TEXT)
UNIQUE_DATA_FILES:=$(shell echo $(call config_column,test_case.config,2) | tr " " "\n" | uniq | tr "\n" " ")
TEMP_DATA_FILE:=build/_tmp.raw

PATTERN_IDS:=$(call config_ids,patterns.config)
UNIQUE_SP_LENGTHS:=$(shell echo $(call config_column,patterns.config,2) | tr " " "\n" | uniq | tr "\n" " ")

ALGOS:=$(call config_ids,algorithms.config)

GM_WORDS = $(foreach TC_ID,$(TC_IDS),$(foreach SP_LENGTH,$(UNIQUE_SP_LENGTHS),collections/$(TC_ID)/patterns/words.$(TC_ID).$(SP_LENGTH).words.txt))
GM_PATTERNS = $(foreach TC_ID,$(TC_IDS),$(foreach PATTERN_ID,$(PATTERN_IDS),collections/$(TC_ID)/patterns/regex.$(TC_ID).$(PATTERN_ID).regex.txt))
GM_EXECS = $(foreach ALGO,$(ALGOS),build/gm_index-$(ALGO).x) \
           $(foreach ALGO,$(ALGOS),build/gm_search-$(ALGO).x) \
           build/create_collection.x
INDEX_SENTINELS = $(foreach TC_ID,$(TC_IDS),$(foreach ALGO,$(ALGOS),collections/$(TC_ID).$(ALGO).sentinel))

RES_FILES = $(foreach TC_ID,$(TC_IDS),$(foreach PATTERN_ID,$(PATTERN_IDS),$(foreach ALGO,$(ALGOS),results/$(ALGO).$(PATTERN_ID).$(TC_ID))))
RES_FILE=results/all.txt

execs: $(GM_EXECS)

collections/%/text.TEXT: build/create_collection.x $(UNIQUE_DATA_FILES)
	$(eval TC_ID:=$(call dim,1,$*))
	$(eval SOURCE_FILE:=$(call config_select,test_case.config,$(TC_ID),2))
	$(eval NUM_BYTES:=$(call config_select,test_case.config,$(TC_ID),5))
	$(eval PREFIX_BYTES:=$(call config_select,test_case.config,$(TC_ID),6))
	@echo "Creating collection $(TC_ID) from file $(SOURCE_FILE)"
	$(if $(findstring -,$(PREFIX_BYTES)), , @head -c $(PREFIX_BYTES) $(SOURCE_FILE) > $(TEMP_DATA_FILE))
	$(if $(findstring -,$(PREFIX_BYTES)), , $(eval SOURCE_FILE:=$(TEMP_DATA_FILE)))
	@build/create_collection.x -i $(SOURCE_FILE) -n $(NUM_BYTES) -c collections/$(TC_ID)
	$(if $(findstring -,$(PREFIX_BYTES)), , @rm $(TEMP_DATA_FILE))

timing: $(RES_FILES)
	cat $(RES_FILES) > $(RES_FILE)
	@cd visualize; make	

# Format: results/[ALGO].[PATTERN_ID].[TC_ID]
results/%: $(GM_EXECS) $(TC_COLLECTIONS) $(GM_PATTERNS)
	@mkdir -p results
	$(eval ALGO:=$(call dim,1,$*))
	$(eval PATTERN_ID:=$(call dim,2,$*))
	$(eval TC_ID:=$(call dim,3,$*))
	$(eval IS_TEXT:=$(call config_select,test_case.config,$(TC_ID),7))
	@echo "# COLL_ID = $(TC_ID)" > $@
	@echo "# PATT_SAMPLE = $(PATTERN_ID)" >> $@
	@echo "# ALGO = $(ALGO)" >> $@
	$(eval SENTINEL:="collections/$(TC_ID).$(ALGO).sentinel")
	@test -e $(SENTINEL) || { echo "Running build/gm_index-$(ALGO).x -c collections/$(TC_ID)"; ./build/gm_index-$(ALGO).x -c collections/$(TC_ID); touch $(SENTINEL); }
	$(MEASURE_CMD) build/gm_search-$(ALGO).x -c collections/$(TC_ID) -p collections/$(TC_ID)/patterns/regex.$(TC_ID).$(PATTERN_ID).regex.txt -t $(IS_TEXT) >> $@
	mv mem-mon-out.csv $@.mem.csv
	Rscript memvis.R $@.mem.csv

../../examples/generate-pattern.x:
	@cd ../../examples; make generate-pattern.x

# Format: collections/[TC_ID]/patterns/words.[TC_ID].[SP_LENGTH].words.txt
collections/%.words.txt: $(TC_COLLECTIONS) ../../examples/generate-pattern.x
	$(eval TC_ID:=$(call dim,2,$*))
	$(eval SP_LENGTH:=$(call dim,3,$*))
	$(eval IS_TEXT:=$(call config_select,test_case.config,$(TC_ID),7))
	@echo "Extracting words of length $(SP_LENGTH) from $(TC_ID)..."
	@../../examples/generate-pattern.x collections/$(TC_ID)/text.TEXT 20 $(SP_LENGTH) $(IS_TEXT) > $@

# Format: collections/[TC_ID]/patterns/regex.[TC_ID].[PATTERN_ID].regex.txt
collections/%.regex.txt: $(GM_WORDS)
	$(eval TC_ID:=$(call dim,2,$*))
	$(eval PATTERN_ID:=$(call dim,3,$*))
	$(eval SP_LENGTH:=$(call config_select,patterns.config,$(PATTERN_ID),2))
	$(eval GAP:=$(call config_select,patterns.config,$(PATTERN_ID),3))
	@echo "$(PATTERN_ID)"
	@echo "Creating regex patterns with words of length $(SP_LENGTH) from $(TC_ID) and a gap of $(GAP)..."
	@sed 'N;s/\n/.{$(GAP)}/' collections/$(TC_ID)/patterns/words.$(TC_ID).$(SP_LENGTH).words.txt > $@

build/%:
	@mkdir -p build
	cd build && cmake -DMEASURE_ENERGY=$(MEASURE_ENERGY) .. && make -j 2

clean-build:
	@echo "Remove executables"
	rm -rf build/*

#clean:
#	rm -f $(GM_PATTERNS)
#	rm -f $(GM_WORDS)
#	rm -f $(INDEX_SENTINELS)

clean-results: 
	rm -f $(RES_FILES) 
		

#cleanall: clean clean-results clean-build
