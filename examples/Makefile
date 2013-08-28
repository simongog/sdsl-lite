include ../Make.helper
CXX_FLAGS=$(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(INC_DIR) -L$(LIB_DIR) 
CCLIB=-lsdsl -ldivsufsort -ldivsufsort64 
SOURCES=$(wildcard *.cpp)
EXECS=$(SOURCES:.cpp=.x)

all: $(EXECS)

build-test: $(EXECS)
	        
%.x:%.cpp 
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(CCLIB) 

clean:
	rm -f $(EXECS)
	rm -rf *.dSYM
