#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <malloc.h>


#include "dacs_basic.h"
#include "stdbool.h"
#include "bitrankw32int.h"

typedef struct sFTRep {
	  uint listLength;
	  t_byte nLevels;
	  uint tamCode;
	  uint * levels;
	  uint * levelsIndex;
	  uint * iniLevel;
	  uint * rankLevels;
	  bitRankW32Int * bS;	
	  //uint * bits_bitmap;
	  uint * base;
	  ushort * base_bits;
	  uint * tablebase;
	  uint tamtablebase;

  	
} FTRep;



// public:
	FTRep* createFT(uint *codewords,uint listLength, uint maxInt);
	uint accessFT(FTRep * listRep,uint param);
	void saveFT(FTRep * listRep, FILE * flist);
	uint * decompressFT(FTRep * listRep, uint n);
	FTRep* loadFT(FILE * flist);
	void destroyFT(FTRep * listRep);
	bool equalsFT(FTRep *lhs, FTRep *rhs);
