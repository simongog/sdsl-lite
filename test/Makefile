include ../Make.helper
CXX_FLAGS= $(MY_CXX_FLAGS) -Wall -Werror -Wunused-parameter -g -O3 -I$(INC_DIR) -L$(LIB_DIR) 
CCLIB=-lsdsl -ldivsufsort -ldivsufsort64 -lgtest 
TMP_DIR=tmp
SOURCES=$(wildcard *Test.cpp)
EXECS=$(SOURCES:.cpp=.x)
EXEC_LIST=$(patsubst %,./%;,$(EXECS))                # list of executables
EXEC_LIST_VAL=$(patsubst %,valgrind ./%;,$(EXECS))   # list of executables preceded by valgrind command

BIT_VECTOR_TC_PATH:=$(call config_column,BitVectorTest.config,2)
BIT_VECTOR_TESTS:=$(patsubst %,bit-vector-test/%,$(BIT_VECTOR_TC_PATH))

RANK_SUPPORT_TESTS:=$(patsubst %,rank-support-test/%,$(BIT_VECTOR_TC_PATH))

SELECT_SUPPORT_TESTS:=$(patsubst %,select-support-test/%,$(BIT_VECTOR_TC_PATH))
SELECT_SUPPORT0_TESTS:=$(patsubst %,select-support0-test/%,$(BIT_VECTOR_TC_PATH))

WT_BYTE_TC_PATHS:=$(call config_column,WtByteTest.config,2)
WT_BYTE_TESTS:=$(patsubst %,wt-byte-test/%,$(WT_BYTE_TC_PATHS))

WT_BYTE_TC_PATHS:=$(call config_column,WtByteTest.config,2)
WT_BYTE_TESTS:=$(patsubst %,wt-byte-test/%,$(WT_BYTE_TC_PATHS))

WT_INT_TC_PATHS:=$(call config_column,WtIntTest.config,2)
WT_INT_TESTS:=$(patsubst %,wt-int-test/%,$(WT_INT_TC_PATHS))

CSA_BYTE_TC_PATHS:=$(call config_column,CsaByteTest.config,2)
CSA_BYTE_TESTS:=$(patsubst %,csa-byte-test/%,$(CSA_BYTE_TC_PATHS))

CSA_INT_TC_PATHS:=$(call config_column,CsaByteTest.config,2)
CSA_INT_TESTS:=$(patsubst %,csa-int-test/%,$(CSA_INT_TC_PATHS))

CST_BYTE_TC_PATHS:=$(call config_column,CstByteTest.config,2)
CST_BYTE_TESTS:=$(patsubst %,cst-byte-test/%,$(CST_BYTE_TC_PATHS))

CST_INT_TC_PATHS:=$(call config_column,CstIntTest.config,2)
CST_INT_TESTS:=$(patsubst %,cst-int-test/%,$(CST_INT_TC_PATHS))

RMQ_TC_PATHS:=$(call config_column,RMQTest.config,2)
RMQ_TESTS:=$(patsubst %,rmq-test/%,$(RMQ_TC_PATHS))

LCP_CONSTRUCT_PATHS:=$(call config_column,LcpConstructTest.config,2)
LCP_CONSTRUCT_TESTS:=$(patsubst %,lcp-construct-test/%,$(LCP_CONSTRUCT_PATHS))

test: bits-test \
      int-vector-test \
      int-vector-buffer-test \
	  bit-vector-test \
	  rank-support-test\
	  select-support-test\
	  select-support0-test\
      wt-byte-test \
      wt-int-test \
      csa-byte-test \
      csa-int-test \
      cst-byte-test \
      cst-int-test \
      rmq-test 

build-test: $(EXECS) 
	cd ../tutorial; make build-test
	cd ../examples; make build-test

build-test-clean: clean
	cd ../tutorial; make clean
	cd ../examples; make clean

generators: BitVectorGenerator.x IntVectorGenerator.x

bits-test: ./BitsTest.x
	./BitsTest.x

int-vector-test: ./IntVectorTest.x
	./IntVectorTest.x

int-vector-buffer-test: ./IntVectorBufferTest.x
	./IntVectorBufferTest.x

bit-vector-test: generators $(BIT_VECTOR_TESTS)

rank-support-test: generators $(RANK_SUPPORT_TESTS)

select-support-test: generators $(SELECT_SUPPORT_TESTS)

select-support0-test: generators $(SELECT_SUPPORT_TESTS)

wt-byte-test: generators $(WT_BYTE_TESTS)

wt-int-test: generators $(WT_INT_TESTS)

csa-byte-test: generators $(CSA_BYTE_TESTS)

csa-int-test: generators $(CSA_INT_TESTS)

cst-byte-test: generators $(CST_BYTE_TESTS)

cst-int-test: generators $(CST_INT_TESTS)

rmq-test: generators $(RMQ_TESTS)

lcp-construct-test: generators $(LCP_CONSTRUCT_TESTS)

bit-vector-test/test_cases/%: ./BitVectorTest.x test_cases/%
	@echo "TEST_CASE: test_cases/$*"
	@$(PREFIX) ./BitVectorTest.x test_cases/$*

rank-support-test/test_cases/%: ./RankSupportTest.x test_cases/%
	@echo "TEST_CASE: test_cases/$*"
	@$(PREFIX) ./RankSupportTest.x test_cases/$*

select-support-test/test_cases/%: ./SelectSupportTest.x test_cases/%
	@echo "TEST_CASE: test_cases/$*"
	@$(PREFIX) ./SelectSupportTest.x test_cases/$*

select-support0-test/test_cases/%: ./SelectSupportTest0.x test_cases/%
	@echo "TEST_CASE: test_cases/$*"
	@$(PREFIX) ./SelectSupport0Test.x test_cases/$*

wt-byte-test/test_cases/%: ./WtByteTest.x test_cases/%
	@echo "TEST_CASE: test_cases/$*"
	@$(PREFIX) ./WtByteTest.x test_cases/$* $(TMP_DIR)/wt-byte.tmp

wt-int-test/test_cases/%: ./WtIntTest.x test_cases/%
	@echo "TEST_CASE: test_cases/$*"
	@$(PREFIX) ./WtIntTest.x test_cases/$* $(TMP_DIR)/wt-int.tmp

csa-byte-test/test_cases/%: ./CsaByteTest.x test_cases/%
	@echo "TEST_CASE: test_cases/$*"
	@$(PREFIX) ./CsaByteTest.x test_cases/$* $(TMP_DIR)/csa-byte.tmp $(TMP_DIR)

csa-int-test/test_cases/%: ./CsaIntTest.x test_cases/%
	$(eval NUM_BYTES:=$(call config_filter,CsaIntTest.config,test_cases/$*,4))
	@echo "TEST_CASE: test_cases/$* NUM_BYTES=$(NUM_BYTES)"
	@$(PREFIX) ./CsaIntTest.x test_cases/$* $(NUM_BYTES) $(TMP_DIR)/csa-int.tmp $(TMP_DIR)

cst-byte-test/test_cases/%: ./CstByteTest.x test_cases/%
	@echo "TEST_CASE: test_cases/$*"
	@$(PREFIX) ./CstByteTest.x test_cases/$* $(TMP_DIR)/cst-byte.tmp $(TMP_DIR)

cst-int-test/test_cases/%: ./CstIntTest.x test_cases/%
	$(eval NUM_BYTES:=$(call config_filter,CstIntTest.config,test_cases/$*,4))
	@echo "TEST_CASE: test_cases/$* NUM_BYTES=$(NUM_BYTES)"
	@$(PREFIX) ./CstIntTest.x test_cases/$* $(NUM_BYTES) $(TMP_DIR)/cst-int.tmp $(TMP_DIR)

rmq-test/test_cases/%: ./RMQTest.x test_cases/%
	@echo "TEST_CASE: test_cases/$*"
	@$(PREFIX) ./RMQTest.x test_cases/$* $(TMP_DIR)/rmq.tmp

lcp-construct-test/test_cases/%: ./LcpConstructTest.x test_cases/%
	$(eval TC_ID:=$(call config_filter,LcpConstructTest.config,test_cases/$*,1))
	@echo "TEST_CASE: test_cases/$*; TC_ID=$(TC_ID)"
	@$(PREFIX) ./LcpConstructTest.x test_cases/$* $(TMP_DIR) $(TC_ID)

# Format: test_cases/int-vec.[LEN].[WIDTH].[DEFAULT]
test_cases/int-vec.%: IntVectorGenerator.x
	@echo "Generate input test_cases/int-vec.$*"
	$(eval LEN:=$(call dim,1,$*)) 
	$(eval WIDTH:=$(call dim,2,$*)) 
	$(eval DEFAULT:=$(call dim,3,$*)) 
	$(eval SEED:=$(call dim,4,$*)) 
	@./IntVectorGenerator.x $@ $(LEN) $(WIDTH) $(DEFAULT) $(SEED)

test_cases/bit-vec.%: BitVectorGenerator.x
	@echo "Generate input test_cases/bit-vec.$*"
	@./BitVectorGenerator.x $@ $*

%.x:%.cpp 
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
	rm -f $(EXECS)
	rm -rf *.dSYM
#	rm -f $(WT_INT_TC_PATHS)

