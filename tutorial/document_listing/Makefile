include ../../Make.helper
CXX_FLAGS=$(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(INC_DIR) -L$(LIB_DIR) 
CCLIB=-lsdsl -ldivsufsort -ldivsufsort64 
EXECS=doc_list_index_sada.x doc_list_index_greedy.x doc_list_index_qprobing.x

all: $(EXECS)

build-test: $(EXECS)

doc_list_index_sada.x: doc_list_index.cpp doc_list_index_sada.hpp
	$(MY_CXX) \
		$(CXX_FLAGS)\
		-DIDX_TYPE="doc_list_index_sada<>" \
		-DIDX_SUF=".sada-idx" \
	-o $@ $< $(CCLIB) 

doc_list_index_greedy.x: doc_list_index.cpp doc_list_index_sada.hpp
	$(MY_CXX) \
		$(CXX_FLAGS)\
		-DIDX_TYPE="doc_list_index_greedy<>" \
		-DIDX_SUF=".greedy-idx" \
	-o $@ $< $(CCLIB) 

doc_list_index_qprobing.x: doc_list_index.cpp doc_list_index_sada.hpp
	$(MY_CXX) \
		$(CXX_FLAGS)\
		-DIDX_TYPE="doc_list_index_qprobing<>" \
		-DIDX_SUF=".greedy-idx" \
	-o $@ $< $(CCLIB) 

clean:
	rm -f $(EXECS)
	rm -rf *.dSYM

