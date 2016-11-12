#ifndef BitRankW32Int_h
#define BitRankW32Int_h
#include "basic.h"
/* bitarray.h
   Copyright (C) 2005, Rodrigo Gonzalez, all rights reserved.

   New RANK, SELECT, SELECT-NEXT and SPARSE RANK implementations.

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/



/////////////
//Rank(B,i)//
/////////////
//_factor = 0  => s=W*lgn
//_factor = P  => s=W*P
//Is interesting to notice
//factor=2 => overhead 50%
//factor=3 => overhead 33%
//factor=4 => overhead 25%
//factor=20=> overhead 5%

typedef struct sbitRankW32Int{
    uint *data;
    char owner;
    uint integers;
    uint factor,b,s;
    uint *Rs;  					//superblock array
    uint n;                  
} bitRankW32Int;
                                 //uso interno para contruir el indice rank
    uint buildRankSub(bitRankW32Int * br,uint ini,uint fin);
    void buildRank(bitRankW32Int * br);            //crea indice para rank

    bitRankW32Int * createBitRankW32Int(uint *bitarray, uint n, char owner, uint factor);
    void destroyBitRankW32Int(bitRankW32Int * br);            //destructor
    uint isBitSet(bitRankW32Int * br, uint i);
    uint rank(bitRankW32Int * br, uint i);           //Nivel 1 bin, nivel 2 sec-pop y nivel 3 sec-bit
    uint lenght_in_bits(bitRankW32Int * br);
    uint prev(bitRankW32Int * br, uint start);       // gives the largest index i<=start such that IsBitSet(i)=true
    uint bselect(bitRankW32Int * br, uint x);         // gives the position of the x:th 1.
    uint select0(bitRankW32Int * br, uint x);        // gives the position of the x:th 0.
    uint select1(bitRankW32Int * br, uint x);        // gives the position of the x:th 1.
    uint spaceRequirementInBits(bitRankW32Int * br);
    /*load-save functions*/
    int save(bitRankW32Int * br, FILE *f);
    int load(bitRankW32Int * br, FILE *f);
    bitRankW32Int * createBitRankW32IntFile(FILE *f, int *error);


#endif
