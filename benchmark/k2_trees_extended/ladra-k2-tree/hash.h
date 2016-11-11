/* Searches using  Horspool's algorithm adapted to 
search inside a text compressed with End-Tagged Dense Code.
Lightweight Natural Language Text Compression: Information Retrieval 2006

Programmed by Antonio Fariña.

Author's contact: Antonio Fariña, Databases Lab, University of
A Coruña. Campus de Elviña s/n. Spain  fari@udc.es

Copyright (C) 2006  Nieves R. Brisaboa, Gonzalo Navarro, Antonio Fariña and José R. Paramá
Author's contact: antonio.fari@gmail.com

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
aint with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/


/*-----------------------------------------------------------------------
 Hash.h: Definition of a HashTable that uses linear hash
 ------------------------------------------------------------------------*/

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <malloc.h>

#include "MemoryManager.h"

#define JUMP 101  		 //jump done when a collision appears
#define OCUP_HASH 1.5	 	 //index of occupation of the hash table
#define SMALL_PRIME 1009 // a small prime number, used to compute a hash function
#define SEED	1159241
/* Type definitions */

#define MIN(a,b) (a < b) ? a : b


struct Nword {
	  unsigned char *word;
	  unsigned int size;
	  unsigned int len;
	  unsigned long weight;
	  unsigned long codeword;
};

typedef struct Nword t_word;

// private:

	MemoryManager _memMgr; 	  /* Holds dynamic memory reserve for words. */

	unsigned long TAM_HASH;   /* # entries in the hash table    */
	unsigned long NumElem;    /* # elements already added to the hash table*/

	unsigned long initialize_hash (unsigned long tamArq);
	unsigned long NearestPrime(unsigned long n);
	unsigned long hashFunction (const unsigned char *aWord, unsigned int len);

// public:
	t_word  *hash;     		  /* holds a hashTable of words  */
	
	unsigned long initialize_hash  (unsigned long tamArq);
	void freeHashTable();
	unsigned long insertElement (const unsigned char *aWord, register unsigned int len,
										register unsigned long *addr);
	unsigned long search (const unsigned char *aWord, register unsigned len,
									unsigned long *returnedAddr);


