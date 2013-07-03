include ../Make.helper
CC=g++
CC_FLAGS=-std=c++11 -Wall -g -O3 -I$(INC_DIR) -L$(LIB_DIR) -DNDEBUG
CCLIB=-lsdsl -ldivsufsort -ldivsufsort64 
SOURCES=$(wildcard *.cpp)
EXECS=$(SOURCES:.cpp=.x)

all: $(EXECS)
	        
%.x:%.cpp
	$(CC) $(CC_FLAGS) -o $@ $< $(CCLIB) 

clean:
	rm -f $(EXECS)
	rm -rf *.dSYM

