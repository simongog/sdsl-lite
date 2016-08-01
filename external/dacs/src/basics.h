
 // Basics
 
#ifndef BASICSINCLUDED
#define BASICSINCLUDED

  // Includes 

#include <sys/types.h>
#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/times.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>



  // Data types

#ifndef byte
	#define byte unsigned char
#endif 

//typedef unsigned char byte;
// typedef unsigned int uint;

//typedef int bool;
//#define true 1
//#define false 0

#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))
	
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define LINE_UP(x,numBytes) ((x>>(4-numBytes)*8))   //move "4-bytes" bytes to the right
						   						   	//so we can write "bytes" bytes from &x
#ifdef LITTLE_ENDIAN
	#define LE_LONG32(x)  (((x&0xff000000)>>24)|((x&0x00ff0000)>>8)|((x&0x0000ff00)<<8)|((x&0x000000ff)<<24))
#else
	#define LE_LONG32(x) x
#endif	

#define mask31 0x0000001F


/* W-1 */
#define Wminusone 31
/*numero de bits del entero de la maquina*/
#define WW 64
/*bits para hacer la mascara para contar mas rapido*/
#define bitsM 8
/*bytes que hacen una palabra */
#define BW 4

  // Bitstream management

//#define W (8*sizeof(uint))
#define W (32)

	// bits needed to represent a number between 0 and n
uint bits (uint n);
        // returns e[p..p+len-1], assuming len <= W
uint bitread (uint *e, uint p, uint len);
        // writes e[p..p+len-1] = s, assuming len <= W
void bitwrite (uint *e, uint p, uint len, uint s);
    // writes e[p..p+len-1] = 0, no assumption on len
        
    /**/ //FARI. WITH ASSUMPTION ON LEN, OR IT CRASHES 
         //NOt WORKING UPON THE LIMIT OF THE STARTING uint.
void bitzero (uint *e, uint p, uint len);
	// reads bit p from e
#define bitget(e,p) (((e)[(p)/W] >> ((p)%W)) & 1)
	// sets bit p in e
#define bitset(e,p) ((e)[(p)/W] |= (1<<((p)%W)))
	// cleans bit p in e
#define bitclean(e,p) ((e)[(p)/W] &= ~(1<<((p)%W)))



/* bitRead and bitWrite as MACROS */
	// returns e[p..p+len-1], assuming len <= W
	//mybitread (uint returned value, uint *e, uint p, uint len)
#define mybitread(answ, v, p, len) \
   { uint *e ; \
   	 e=v;\
     e += p/W; p %= W; \
     answ = *e >> p; \
     if (len == W) \
	  { if (p) answ |= (*(e+1)) << (W-p); \
	  } \
     else { if (p+len > W) answ |= (*(e+1)) << (W-p); \
            answ &= (1<<len)-1; \
	  } \
   }


  	// writes e[p..p+len-1] = s, len <= W
	//void bitwrite (uint *e, uint p, uint len, uint s)
#define mybitwrite(v, p, len, s) \
   { uint *e ; \
   	 e=v; \
    e += p/W; p %= W; \
     if (len == W) \
	  { *e |= (*e & ((1<<p)-1)) | (s << p); \
            if (p) { \
            	e++; \
            	*e = (*e & ~((1<<p)-1)) | (s >> (W-p)); \
        	} \
	  } \
     else { if (p+len <= W) \
	       { *e = (*e & ~(((1<<len)-1)<<p)) | (s << p); \
	       } \
	       else { \
	    	*e = (*e & ((1<<p)-1)) | (s << p); \
            e++; len -= W-p; \
            *e = (*e & ~((1<<len)-1)) | (s >> (W-p)); \
           } \
	  } \
   } 

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
