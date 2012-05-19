GCCFLAGS=-g -O3 -DNDEBUG -funroll-loops
GCCINC=-I${HOME}/include/ 
GCCLIBDIR=-L${HOME}/lib/ 
GCCLIB=-lsdsl -ldivsufsort -ldivsufsort64 -lgtest -lpthread -msse4.2


all: IntVectorTest RankSupportTest SelectSupportTest\
	 SelectSupportTest0 BitMagicTest BitMagicTestSSE \
	 BitVectorTest \
	 CstTest

IntVectorTest: IntVectorTest.cpp
	g++ ${GCCFLAGS} ${GCCINC} ${GCCLIBDIR} IntVectorTest.cpp -o IntVectorTest ${GCCLIB}

RankSupportTest: RankSupportTest.cpp
	g++ ${GCCFLAGS} ${GCCINC} ${GCCLIBDIR} RankSupportTest.cpp -o RankSupportTest ${GCCLIB}

SelectSupportTest: SelectSupportTest.cpp
	g++ ${GCCFLAGS} ${GCCINC} ${GCCLIBDIR} SelectSupportTest.cpp -o SelectSupportTest ${GCCLIB}

SelectSupportTest0: SelectSupportTest0.cpp
	g++ ${GCCFLAGS} ${GCCINC} ${GCCLIBDIR} SelectSupportTest0.cpp -o SelectSupportTest0 ${GCCLIB}

BitMagicTest: BitMagicTest.cpp
	g++ ${GCCFLAGS} ${GCCINC} ${GCCLIBDIR} BitMagicTest.cpp -o BitMagicTest ${GCCLIB}

BitMagicTestSSE: BitMagicTest.cpp
	g++ -msse4.2 ${GCCFLAGS} ${GCCINC} ${GCCLIBDIR} BitMagicTest.cpp -o BitMagicTestSSE ${GCCLIB}

BitVectorTest: BitVectorTest.cpp
	g++ -msse4.2 ${GCCFLAGS} ${GCCINC} ${GCCLIBDIR} BitVectorTest.cpp -o BitVectorTest ${GCCLIB}

CstTest: CstTest.cpp
	g++ ${GCCFLAGS} ${GCCINC} ${GCCLIBDIR} CstTest.cpp -o CstTest ${GCCLIB}

test:
	./IntVectorTest
	./BitVectorTest
	./RankSupportTest
	./SelectSupportTest
	./SelectSupportTest0
	./BitMagicTest
	./BitMagicTestSSE
	./test_cases/small/get_corpus.sh
	./CstTest

clean:
	rm -f IntVectorTest RankSupportTest SelectSupportTest \
		  SelectSupportTest0 BitMagicTest BitMagicTestSSE BitVectorTest \
		  CstTest 
	rm -rf *.dSYM
