include ../Make.helper
CC=g++
CC_FLAGS=-Wall -Werror -Wunused-parameter -g -O3 -I$(INC_DIR) -L$(LIB_DIR) 
CCLIB=-lsdsl -ldivsufsort -ldivsufsort64 -lgtest 
TMP_DIR=tmp
SOURCES=$(wildcard *Test.cpp)
EXECS=$(SOURCES:.cpp=.x)
EXEC_LIST=$(patsubst %,./%;,$(EXECS))                # list of executables
EXEC_LIST_VAL=$(patsubst %,valgrind ./%;,$(EXECS))   # list of executables preceded by valgrind command


WT_INT_TC_PATHS:=$(call config_column,WtIntTest.config,2)
WT_INT_TESTS:=$(patsubst %,%.wt-int-test ,$(WT_INT_TC_PATHS))

WT_BYTE_TC_PATHS:=$(call config_column,WtByteTest.config,2)
WT_BYTE_TESTS:=$(patsubst %,%.wt-byte-test ,$(WT_BYTE_TC_PATHS))

all: $(WT_INT_TC_PATHS) \
     $(WT_BYTE_TC_PATHS) \
     $(EXECS)

wt-byte-test: $(WT_BYTE_TESTS)

wt-int-test: $(WT_INT_TESTS)

test_cases/%.wt-byte-test: ./WtByteTest.x test_cases/%
	@echo "test_cases/$*"
	@./WtByteTest.x test_cases/$* $(TMP_DIR)/wt-byte.tmp

test_cases/%.wt-int-test: ./WtIntTest.x test_cases/%
	./WtIntTest.x test_cases/$* $(TMP_DIR)/wt-int.tmp

%.x:%.cpp /Users/sgog/lib/libsdsl.a
	$(CC) $(CC_FLAGS) -o $@ $< $(CCLIB) 

test_cases/wt-int.%: WtIntTestGenerateTest.x
	$(eval LEN:=$(call dim,1,$*)) 
	$(eval WIDTH:=$(call dim,2,$*)) 
	$(eval DEFAULT:=$(call dim,3,$*)) 
	@./WtIntTestGenerateTest.x $@ $(LEN) $(WIDTH) $(DEFAULT)


test:
	./test_cases/small/get_corpus.sh
	$(EXEC_LIST)

vtest:
	./test_cases/small/get_corpus.sh
	$(EXEC_LIST_VAL)

test_cases/%: 
	$(eval URL:=$(call config_filter,download.config,$@,2))
	@$(if $(URL),,\
		$(error "No downlaod link specified for test case $@") )
	@echo "Downlaod input from $(URL) using curl"
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

