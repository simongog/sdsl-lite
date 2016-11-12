/* Searches using  Horspool's algorithm adapted to 
search inside a text compressed with End-Tagged Dense Code.
Lightweight Natural Language Text Compression: Information Retrieval 2006

Programmed by Antonio Faria.

Author's contact: Antonio Faria, Databases Lab, University of
A Corua. Campus de Elvia s/n. Spain  fari@udc.es

Copyright (C) 2006  Nieves R. Brisaboa, Gonzalo Navarro, Antonio Faria and Jos R. Param
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
 Hash.c: Definition of a HashTable that uses linear hash
 ------------------------------------------------------------------------*/

#include "hash.h"

#define mystrncpy(a,b,n) {unsigned int i; for(i=0;i<n;i++) a[i]=b[i];}

/*-----------------------------------------------------------------
 Initilization of data structures used by the hashTable
 ---------------------------------------------------------------- */
unsigned long initialize_hash (unsigned long sizeVoc) {

	unsigned long i;

	NumElem = 0;
	TAM_HASH = OCUP_HASH * sizeVoc;
	TAM_HASH = NearestPrime(TAM_HASH);
	hash = (t_word *) malloc(TAM_HASH * sizeof(t_word));

	for(i = 0; i < TAM_HASH; i++)    {
		hash[i].word = NULL;
		hash[i].weight = 0;
		hash[i].size = 0;
		hash[i].codeword = 0;
		hash[i].len=0;
	}

	return TAM_HASH;
}

/*-----------------------------------------------------------------
 Frees the memory allocated by the hash table 
-----------------------------------------------------------------*/
void freeHashTable() {
	destroyMemoryManager(_memMgr);
	free(hash);
}

/*------------------------------------------------------------------
 Find the nearest prime number over n. 
 ---------------------------------------------------------------- */
unsigned long NearestPrime(unsigned long n)
{
    long position;  /* the prime number being sought */
    long index;  

    for (position = n; ; position++)
    {
        // checks if those values from 2 to $\sqrt{m}$ can be factors of $m$ */
        for (index = 2; index <= (long) sqrt((double) position) && position % index != 0; index++) ;

        if (position % index != 0)  /* No factors in that range, therefore a prime number was found */
        {
            break;
        }
    }
    return position;
} 

/*------------------------------------------------------------------*
//   http://goanna.cs.rmit.edu.au/~hugh/software/zwh-ipl/
//	/* Author J. Zobel, April 2001.
//	   Permission to use this code is freely granted, provided that this
//	   statement is retained. */
//
//	/* Bitwise hash function.  Note that tsize does not have to be prime. */
//	unsigned int bitwisehash(char *word, int tsize, unsigned int seed)
//	{
//	    char	c;
//	    unsigned int h;
//
//	    h = seed;
//	    for( ; ( c=*word )!='\0' ; word++ )
//	    {
//		h ^= ( (h << 5) + c + (h >> 2) );
//	    }
//	    return((unsigned int)((h&0x7fffffff) % tsize));
//	}
/*---------------------------------------------------------------- */

/*------------------------------------------------------------------
 Modification of Zobel's bitwise function to have into account the 
 lenght of the key explicetely 
 ---------------------------------------------------------------- */
unsigned long hashFunction (const unsigned char *aWord, unsigned int len)
{
    char c;
    register unsigned int h;
    register unsigned long last;
    last=((unsigned long) aWord) +len -1;

    h = SEED;
   
    for( ; ((unsigned long) aWord <=last ) ; )
    {
    	c=*(aWord++);
		//c=*aWord;
		h ^= ( (h << 5) + c + (h >> 2) );
    }
    return((unsigned int)((h&0x7fffffff) % TAM_HASH));
}


/*-----------------------------------------------------------------------
  compares two strings
 ---------------------------------------------------------------------*/

 /*------------------------------------------------------------------*
//   http://goanna.cs.rmit.edu.au/~hugh/software/zwh-ipl/
//	/* Author J. Zobel, April 2001.
//	   Permission to use this code is freely granted, provided that this
//	   statement is retained. */
//
//	/* String compare function. */
//	int scmp( char *s1, char *s2 )
//	{
//	    while( *s1 != '\0' && *s1 == *s2 )
//	    {
//		s1++;
//		s2++;
//	    }
//	    return( *s1-*s2 );
//	}
/*---------------------------------------------------------------- */


/*------------------------------------------------------------------
 Modification of Zobel's scmp function compare two strings
 ---------------------------------------------------------------- */
inline int strcomp(const unsigned char *s1, const unsigned char *s2, register unsigned int ws1, unsigned int ws2) {
	 
	 if (ws1 !=ws2) {
	 		return -1;
	 }
	 
	 {  register unsigned long iters;	    
	 	 iters=1;
	    while( iters<ws1 && *s1 == *s2 )
	    {
			s1++;
			s2++;
			iters++;
	    }
		
	    return( *s1-*s2 );
	 }
}

/*------------------------------------------------------------------
 Modification of Zobel's scmp function compare two strings
 ---------------------------------------------------------------- */
/* inline int strcomp_vOLD(const unsigned char *s1, const unsigned char *s2, unsigned int ws1, unsigned int ws2) {
	 
	 if (ws1 !=ws2) {
	 		return -1;
	 }
	 
	 {  register unsigned int iters;
	    register unsigned int end;
	    end = MIN(ws1,ws2);
	 	 iters=1;
	    while( iters<end && *s1 == *s2 )
	    {
			s1++;
			s2++;
			iters++;
	    }
	    return( *s1-*s2 );
	 }
} */


/*-----------------------------------------------------------------------
 Inserts a new word in a given position of the hashTable (position previously computed)
 ---------------------------------------------------------------------*/
 unsigned long insertElement (const unsigned char *aWord, register unsigned int len,
                                         register unsigned long *addr) {
	
	if(*addr == TAM_HASH)
	{
		printf("Not enough memory, vocabulary exceeds maximun size !\n");
		exit(1);
	}

    getMemoryBlock(_memMgr,( byte **)&(hash[*addr].word),len+1);
   
	//mystrncpy((char *) hash[*addr].word, (char *)aWord, len);

	{
		unsigned int i; 
		for(i=0;i<len;i++) 
			hash[*addr].word[i]=aWord[i];
		}
	hash[*addr].word[len]='\0';
	//hash[*addr].size = len;
	hash[*addr].weight=0;
	NumElem++;

	return *addr;
}

/*-----------------------------------------------------------------------
 Searches for a word in the hash table and returns its position in the
 vocabulary. It returns the next "available" position in the vocabulary,
 if the word is not in the hash table. That is: a "0-node" position.
 It also returns -using attribute returnedAddr- the position where word
 was found (or where it should go if it was inserted in next "insert call".
 -----------------------------------------------------------------------*/
unsigned long search (const unsigned char *aWord, register unsigned len,
								 unsigned long *returnedAddr){

	register unsigned long addr, Saddr;

	addr = hashFunction(aWord,len);
	
	Saddr = addr;


	while((hash[addr].word  != NULL)&&((strcomp(hash[addr].word, aWord,  hash[addr].len, len)) != 0))  {
		addr = ((addr + JUMP) %TAM_HASH);
	}

	*returnedAddr = addr;  // position returned

	if(hash[addr].word  == NULL) {
		return NumElem;	//Word was not found
	}
	
	return 0; 			//Word was found	
}

