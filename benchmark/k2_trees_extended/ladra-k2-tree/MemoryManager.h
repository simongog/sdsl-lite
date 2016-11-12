/* Searches using text compressed and indexed with WTDC.

Programmed by Susana Ladra.

Author's contact: Susana Ladra, Databases Lab, University of
A Coruna. Campus de Elvina s/n. Spain  sladra@udc.es

Copyright (C) 2006  Nieves R. Brisaboa, Antonio Farina and Susana Ladra
Author's contact: asusanaladra@gmail.com

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/


/*-----------------------------------------------------------------------
MemoryManager.h:     Purpose    : Improve hash performance.
Reserves large blocks of memory and gives pointers to small  portions of that block 
when requested. This improves performance since a unique "LARGE ALLOCATION"
of memory is needed (a unique call to malloc). It is also responsible of freeing memory.  
 ------------------------------------------------------------------------*/

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include <malloc.h>

#ifndef DEFBYTE	
	#define byte unsigned char
	#define DEFBYTE
#endif

#define LARGE_BLOCK_SIZE 1048576 // Size of the blocks of memory that will be allocated
#define MAX_BLOCKS 1024				// Maximum number of blocks of size LARGE_BLOCK_SIZE that
											// can be allocated 

	/*
	 * Definition of structure MemoryManager
	 */
	struct sMem {
		byte *BLOCKS[MAX_BLOCKS]; 		//array of blocks of size LARGE_BLOCK_SIZE
		unsigned int  currentBlock;	//currentBlock in the array of blocks
		unsigned long remainderBytes; //number of bytes not yet assigned in BLOCKS[currentBlock]
		byte *availableByte;				//pointer to next byte not yet assigned 
	} ;
	
	typedef struct sMem *MemoryManager;
	
		 
	MemoryManager createMemoryManager(void);
	void destroyMemoryManager (MemoryManager mm);
	void getMemoryBlock (MemoryManager mm, byte **dst, const unsigned int size);
	void createNewMemoryBlock (MemoryManager mm);

