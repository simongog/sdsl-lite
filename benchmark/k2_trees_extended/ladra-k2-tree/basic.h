#ifndef BASICSINCLUDED
#define BASICSINCLUDED

#include <sys/types.h>
#include <sys/resource.h>
#include <sys/times.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>

/* basics.h
   Copyright (C) 2005, Rodrigo Gonzalez, all rights reserved.

   Some preliminary stuff

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



#define mask31 0x0000001F

#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))

/*numero de bits del entero de la maquina*/
#define W 32
/* W-1 */
#define Wminusone 31
/*numero de bits del entero de la maquina*/
#define WW 64
/*bits para hacer la mascara para contar mas rapido*/
#define bitsM 8
/*bytes que hacen una palabra */
#define BW 4
#ifndef uchar
#define uchar unsigned char
#endif
#ifndef uint
#define uint unsigned int
#endif
#ifndef ulong
#define ulong unsigned long
#endif
#define size_uchar 256

#ifndef byte
      #define byte unsigned char
#endif

/* reads bit p from e */
#define bitgetchar(e,p) ((((e)[(p)/bitsM] >> ((p)%bitsM))) & 1)
/* sets bit p in e */
#define bitsetchar(e,p) ((e)[(p)/bitsM] |= (1<<((p)%bitsM)))
/* cleans bit p in e */
#define bitcleanchar(e,p) ((e)[(p)/bitsM] &= ~(1<<((p)%bitsM)))


/* reads bit p from e */
#define bitget(e,p) ((((e)[(p)/W] >> ((p)%W))) & 1)
/* sets bit p in e */
#define bitset(e,p) ((e)[(p)/W] |= (1<<((p)%W)))
/* cleans bit p in e */
#define bitclean(e,p) ((e)[(p)/W] &= ~(1<<((p)%W)))

/* numero de enteros necesarios para representar e elementos de largo n */
#define enteros(e,n) ((e)*(n))/W+(((e)*(n))%W > 0)
/* bits needed to represent a number between 0 and n */
uint bits (uint n);
uint bitread (uint *e, uint p, uint len);
        // writes e[p..p+len-1] = s, assuming len <= W
void bitwrite (uint *e, uint p, uint len, uint s);
    // writes e[p..p+len-1] = 0, no assumption on len


uint GetField(uint *A, register  uint len, register uint index);
void SetField(uint *A,register uint len, register uint index,register  uint x);


uint GetVarField(uint *A, register  uint ini, register uint fin);

void SetVarField(uint *A,register uint ini, register uint fin,register  uint x);

unsigned GetFieldW32(uint *A,register uint index);


void SetField32(uint *A, register uint index,register uint x) ;

unsigned GetFieldW16(uint *A,register uint index);
unsigned GetFieldW4(uint *A,register uint index) ;

uint popcount (register int x);

uint popcount16 (register int x);

uint popcount8 (register int x);
#endif
