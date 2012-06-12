GCCFLAGS=-g -O3 -DNDEBUG -funroll-loops
GCCINC=-I${HOME}/include/ 
GCCLIBDIR=-L${HOME}/lib/ 
GCCLIB=-lsdsl -ldivsufsort -ldivsufsort64 -msse4.2


all: cst_bfs_iterator louds_tree wt_huff_from_text json2html


cst_bfs_iterator: cst_bfs_iterator.cpp
	g++ ${GCCFLAGS} ${GCCINC} ${GCCLIBDIR} cst_bfs_iterator.cpp -o cst_bfs_iterator ${GCCLIB}

louds_tree: louds_tree.cpp
	g++ ${GCCFLAGS} ${GCCINC} ${GCCLIBDIR} louds_tree.cpp -o louds_tree ${GCCLIB}

wt_huff_from_text: wt_huff_from_text.cpp
	g++ ${GCCFLAGS} ${GCCINC} ${GCCLIBDIR} wt_huff_from_text.cpp -o wt_huff_from_text ${GCCLIB}

json2html: json2html.cpp
	g++ ${GCCFLAGS} ${GCCINC} ${GCCLIBDIR} json2html.cpp -o json2html ${GCCLIB}

clean:
	rm -f cst_bfs_iterator louds_tree wt_huff_from_text json2html
	rm -rf *.dSYM
