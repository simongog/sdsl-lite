include ../Make.helper
CXX_FLAGS=$(MY_CXX_FLAGS) -Wall -Werror -Wunused-parameter -g -O3 \
		  -I$(INC_DIR) \
		  -I../build/external/gtest-1.6.0/include \
		  -L$(LIB_DIR) \
		  -L../build/external/gtest-1.6.0
LIB_SDSL=$(LIB_DIR)/libsdsl.a
CCLIB=-lsdsl -ldivsufsort -ldivsufsort64 -lgtest
TMP_DIR=tmp
SOURCES=$(wildcard *Test.cpp)
EXECS=$(SOURCES:.cpp=.x)
EXEC_LIST=$(patsubst %,./%;,$(EXECS))                # list of executables

GCDAs=$(SOURCES:.cpp=.gcda)
GCNOs=$(SOURCES:.cpp=.gcno)

BIT_VECTOR_TC_PATH:=$(call config_column,BitVectorTest.config,2)
BIT_VECTOR_TESTS:=$(patsubst %,bit-vector-test/%,$(BIT_VECTOR_TC_PATH))

RANK_SUPPORT_TESTS:=$(patsubst %,rank-support-test/%,$(BIT_VECTOR_TC_PATH))

SELECT_SUPPORT_TESTS:=$(patsubst %,select-support-test/%,$(BIT_VECTOR_TC_PATH))

WT_BYTE_TC_PATHS:=$(call config_column,WtByteTest.config,2)
WT_BYTE_TESTS:=$(patsubst %,wt-byte-test/%,$(WT_BYTE_TC_PATHS))

WT_INT_TC_PATHS:=$(call config_column,WtIntTest.config,2)
WT_INT_TESTS:=$(patsubst %,wt-int-test/%,$(WT_INT_TC_PATHS))

CSA_BYTE_TC_PATHS:=$(call config_column,CsaByteTest.config,2)
CSA_BYTE_TESTS:=$(patsubst %,csa-byte-test/%,$(CSA_BYTE_TC_PATHS))

CSA_INT_TC_PATHS:=$(call config_column,CsaIntTest.config,2)
CSA_INT_TESTS:=$(patsubst %,csa-int-test/%,$(CSA_INT_TC_PATHS))

CST_BYTE_TC_PATHS:=$(call config_column,CstByteTest.config,2)
CST_BYTE_TESTS:=$(patsubst %,cst-byte-test/%,$(CST_BYTE_TC_PATHS))

CST_INT_TC_PATHS:=$(call config_column,CstIntTest.config,2)
CST_INT_TESTS:=$(patsubst %,cst-int-test/%,$(CST_INT_TC_PATHS))

RMQ_TC_PATHS:=$(call config_column,RMQTest.config,2)
RMQ_TESTS:=$(patsubst %,rmq-test/%,$(RMQ_TC_PATHS))

LCP_CONSTRUCT_PATHS:=$(call config_column,LcpConstructTest.config,2)
LCP_CONSTRUCT_TESTS:=$(patsubst %,lcp-construct-test/%,$(LCP_CONSTRUCT_PATHS))

SA_CONSTRUCT_PATHS:=$(call config_column,SaConstructTest.config,2)
SA_CONSTRUCT_TESTS:=$(patsubst %,sa-construct-test/%,$(SA_CONSTRUCT_PATHS))

SEARCH_BIDIRECTIONAL_PATHS:=$(call config_column,SearchBidirectionalTest.config,2)
SEARCH_BIDIRECTIONAL_TESTS:=$(patsubst %,search-bidirectional-test/%,$(SEARCH_BIDIRECTIONAL_PATHS))

TC_PATHS:= $(WT_BYTE_TC_PATHS) $(WT_INT_TC_PATHS) \
		   $(CSA_BYTE_TC_PATHS) $(CSA_INT_TC_PATHS) \
		   $(CST_BYTE_TC_PATHS) $(CST_INT_TC_PATHS)



# do not delete the generated/downloaded test_cases
.SECONDARY: $(TC_PATHS)

test:	bits-test \
	coder-test \
	int-vector-test \
	inv-perm-support-test \
	int-vector-mapper-test \
	int-vector-buffer-test \
	bit-vector-test \
	rank-support-test\
	select-support-test\
	wt-byte-test \
	wt-int-test \
	csa-byte-test \
	csa-int-test \
	cst-byte-test \
	cst-int-test \
	rmq-test \
	sa-construct-test \
	search-bidirectional-test \
	sorted-int-stack-test \
	nn-dict-dynamic-test

build-test: $(EXECS)
	cd ../tutorial; make build-test
	cd ../examples; make build-test
	cd ../benchmark/indexing_count; make
	cd ../benchmark/indexing_locate; make
	cd ../benchmark/indexing_extract; make
	cd ../benchmark/rrr_vector; make
	cd ../benchmark/document_retrieval; make
	cd ../benchmark/wavelet_trees; make

build-test-clean: clean
	cd ../tutorial; make clean
	cd ../examples; make clean
	cd ../benchmark/indexing_count; make clean-build
	cd ../benchmark/indexing_locate; make clean-build
	cd ../benchmark/indexing_extract; make clean-build
	cd ../benchmark/rrr_vector; make clean-build
	cd ../benchmark/document_retrieval; make clean-build
	cd ../benchmark/wavelet_trees; make clean-build

generators: BitVectorGenerator.x IntVectorGenerator.x ReplaceIntVectorValue.x

bits-test: ./BitsTest.x
	@$(PREFIX) ./BitsTest.x

coder-test: ./CoderTest.x
	@$(PREFIX) ./CoderTest.x

int-vector-test: ./IntVectorTest.x
	@$(PREFIX) ./IntVectorTest.x

inv-perm-support-test: ./InvPermSupportTest.x
	@$(PREFIX) ./InvPermSupportTest.x

int-vector-buffer-test: ./IntVectorBufferTest.x
	@$(PREFIX) ./IntVectorBufferTest.x

int-vector-mapper-test: ./IntVectorMapperTest.x
	@$(PREFIX) ./IntVectorMapperTest.x

sorted-int-stack-test: ./SortedIntStackTest.x
	@$(PREFIX) ./SortedIntStackTest.x

nn-dict-dynamic-test: ./NNDictDynamicTest.x
	@$(PREFIX) ./NNDictDynamicTest.x

bit-vector-test: generators $(BIT_VECTOR_TESTS)

rank-support-test: generators $(RANK_SUPPORT_TESTS)

select-support-test: generators $(SELECT_SUPPORT_TESTS)

wt-byte-test: generators $(WT_BYTE_TESTS)

wt-int-test: generators $(WT_INT_TESTS)

csa-byte-test: generators $(CSA_BYTE_TESTS)

csa-int-test: generators $(CSA_INT_TESTS)

cst-byte-test: generators $(CST_BYTE_TESTS)

cst-int-test: generators $(CST_INT_TESTS)

rmq-test: generators $(RMQ_TESTS)

lcp-construct-test: generators $(LCP_CONSTRUCT_TESTS)

sa-construct-test: generators $(SA_CONSTRUCT_TESTS)

search-bidirectional-test: generators $(SEARCH_BIDIRECTIONAL_TESTS)

bit-vector-test/test_cases/%: ./BitVectorTest.x test_cases/%
	@echo "TEST_CASE: test_cases/$*"
	@$(PREFIX) ./BitVectorTest.x test_cases/$*

rank-support-test/test_cases/%: ./RankSupportTest.x test_cases/%
	@echo "TEST_CASE: test_cases/$*"
	@$(PREFIX) ./RankSupportTest.x test_cases/$*

select-support-test/test_cases/%: ./SelectSupportTest.x test_cases/%
	@echo "TEST_CASE: test_cases/$*"
	@$(PREFIX) ./SelectSupportTest.x test_cases/$*

wt-byte-test/test_cases/%: ./WtByteTest.x test_cases/%
	@echo "TEST_CASE: test_cases/$* semi-external"
	@$(PREFIX) ./WtByteTest.x test_cases/$* $(TMP_DIR)/wt-byte.tmp
	@echo "TEST_CASE: test_cases/$* in-memory"
	@$(PREFIX) ./WtByteTest.x test_cases/$* $(TMP_DIR)/wt-byte.tmp in-memory

wt-int-test/test_cases/%: ./WtIntTest.x test_cases/%
	@echo "TEST_CASE: test_cases/$* semi-external"
	@$(PREFIX) ./WtIntTest.x test_cases/$* $(TMP_DIR)/wt-int.tmp
	@echo "TEST_CASE: test_cases/$* in-memory"
	@$(PREFIX) ./WtIntTest.x test_cases/$* $(TMP_DIR)/wt-int.tmp in-memory

csa-byte-test/test_cases/%: ./CsaByteTest.x test_cases/%
	@echo "TEST_CASE: test_cases/$* semi-external"
	@$(PREFIX) ./CsaByteTest.x test_cases/$* $(TMP_DIR)/csa-byte.tmp $(TMP_DIR)
	@echo "TEST_CASE: test_cases/$* in-memory"
	@$(PREFIX) ./CsaByteTest.x test_cases/$* $(TMP_DIR)/csa-byte.tmp $(TMP_DIR) in-memory

csa-int-test/test_cases/%: ./CsaIntTest.x test_cases/%
	$(eval NUM_BYTES:=$(call config_filter,CsaIntTest.config,test_cases/$*,4))
	@echo "TEST_CASE: test_cases/$* NUM_BYTES=$(NUM_BYTES) semi-external"
	@$(PREFIX) ./CsaIntTest.x test_cases/$* $(NUM_BYTES) $(TMP_DIR)/csa-int.tmp $(TMP_DIR)
	@echo "TEST_CASE: test_cases/$* NUM_BYTES=$(NUM_BYTES) in-memory"
	@$(PREFIX) ./CsaIntTest.x test_cases/$* $(NUM_BYTES) $(TMP_DIR)/csa-int.tmp $(TMP_DIR) in-memory

cst-byte-test/test_cases/%: ./CstByteTest.x test_cases/%
	@echo "TEST_CASE: test_cases/$* semi-external"
	$(PREFIX) ./CstByteTest.x test_cases/$* $(TMP_DIR)/cst-byte.tmp $(TMP_DIR)
	@echo "TEST_CASE: test_cases/$* in-memory"
	$(PREFIX) ./CstByteTest.x test_cases/$* $(TMP_DIR)/cst-byte.tmp $(TMP_DIR) in-memory

cst-int-test/test_cases/%: ./CstIntTest.x test_cases/%
	$(eval NUM_BYTES:=$(call config_filter,CstIntTest.config,test_cases/$*,4))
	@echo "TEST_CASE: test_cases/$* NUM_BYTES=$(NUM_BYTES) semi-external"
	@$(PREFIX) ./CstIntTest.x test_cases/$* $(NUM_BYTES) $(TMP_DIR)/cst-int.tmp $(TMP_DIR)
	@echo "TEST_CASE: test_cases/$* NUM_BYTES=$(NUM_BYTES) in-memory"
	@$(PREFIX) ./CstIntTest.x test_cases/$* $(NUM_BYTES) $(TMP_DIR)/cst-int.tmp $(TMP_DIR) in-memory

rmq-test/test_cases/%: ./RMQTest.x test_cases/%
	@echo "TEST_CASE: test_cases/$*"
	@$(PREFIX) ./RMQTest.x test_cases/$* $(TMP_DIR)/rmq.tmp

lcp-construct-test/test_cases/%: ./LcpConstructTest.x test_cases/%
	$(eval TC_ID:=$(call config_filter,LcpConstructTest.config,test_cases/$*,1))
	@echo "TEST_CASE: test_cases/$*; TC_ID=$(TC_ID)"
	@$(PREFIX) ./LcpConstructTest.x test_cases/$* $(TMP_DIR) $(TC_ID)

sa-construct-test/test_cases/%: ./SaConstructTest.x test_cases/%
	$(eval TC_ID:=$(call config_filter,SaConstructTest.config,test_cases/$*,1))
	@echo "TEST_CASE: test_cases/$*; TC_ID=$(TC_ID)"
	@$(PREFIX) ./SaConstructTest.x test_cases/$* $(TMP_DIR) $(TC_ID)

search-bidirectional-test/test_cases/%: ./SearchBidirectionalTest.x test_cases/%
	@echo "TEST_CASE: test_cases/$*"
	@$(PREFIX) ./SearchBidirectionalTest.x test_cases/$*

# Format: test_cases/int-vec.[LEN].[WIDTH].[DEFAULT]
test_cases/int-vec.%: IntVectorGenerator.x
	@echo "Generate input test_cases/int-vec.$*"
	$(eval LEN:=$(call dim,1,$*))
	$(eval WIDTH:=$(call dim,2,$*))
	$(eval DEFAULT:=$(call dim,3,$*))
	$(eval SEED:=$(call dim,4,$*))
	@./IntVectorGenerator.x $@ $(LEN) $(WIDTH) $(DEFAULT) $(SEED)

# Format: test_cases/int-vec-sa.[LEN].[WIDTH].[DEFAULT]
test_cases/int-vec-sa.%: IntVectorGenerator.x
	@echo "Generate input test_cases/int-vec-sa.$*"
	$(eval LEN:=$(call dim,1,$*))
	$(eval WIDTH:=$(call dim,2,$*))
	$(eval DEFAULT:=$(call dim,3,$*))
	$(eval SEED:=$(call dim,4,$*))
	@./IntVectorGenerator.x $@ $(LEN) $(WIDTH) $(DEFAULT) $(SEED)
	@./ReplaceIntVectorValue.x $@ 0 1



test_cases/bit-vec.%: BitVectorGenerator.x
	@echo "Generate input test_cases/bit-vec.$*"
	@./BitVectorGenerator.x $@ $*

%.x:%.cpp $(LIB_SDSL)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(CCLIB)

test_cases/%:
	$(eval URL:=$(call config_filter,download.config,$@,2))
	@$(if $(URL),,\
		$(error "No download link specified for test case $@") )
	@echo "Download input from $(URL) using curl"
	$(eval DEST_DIR:=$(shell dirname $@))
	cd $(DEST_DIR); curl -O $(URL)
	$(eval FILE:=$(DEST_DIR)/$(notdir $(URL)))
	@$(if $(filter-out ".gz",$(FILE)),\
		echo "Extract file $(FILE) using gunzip";\
		gunzip $(FILE))


clean:
	rm -f $(EXECS) $(GCDAs) $(GCNOs)
	rm -f tmp/*
	rm -rf *.dSYM

cleanall: clean
	rm -f $(TC_PATHS)
