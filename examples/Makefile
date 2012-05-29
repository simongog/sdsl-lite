GCCFLAGS=-g -O3 -DNDEBUG -funroll-loops
GCCINC=-I${HOME}/include/ 
GCCLIBDIR=-L${HOME}/lib/ 
GCCLIB=-lsdsl -ldivsufsort -ldivsufsort64 -msse4.2


all: cst_bfs_iterator louds_tree wt_huff_from_text


cst_bfs_iterator: cst_bfs_iterator.cpp
	g++ ${GCCFLAGS} ${GCCINC} ${GCCLIBDIR} cst_bfs_iterator.cpp -o cst_bfs_iterator ${GCCLIB}

louds_tree: louds_tree.cpp
	g++ ${GCCFLAGS} ${GCCINC} ${GCCLIBDIR} louds_tree.cpp -o louds_tree ${GCCLIB}

wt_huff_from_text: wt_huff_from_text.cpp
	g++ ${GCCFLAGS} ${GCCINC} ${GCCLIBDIR} wt_huff_from_text.cpp -o wt_huff_from_text ${GCCLIB}



clean:
	rm -f cst_bfs_iterator louds_tree wt_huff_from_text
	rm -rf *.dSYM
