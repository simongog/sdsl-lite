#ifndef KTREE_H
#define KTREE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "directcodes.h"
#include "hash.h"

#include "bitrankw32int.h"
#include "basic.h"
#define MAX_INFO 1024*1024*100+10

#define K1 4
#define K2 2

//#define L 0
//#define K2_2 4
//#define K2_3 16
//#define K2_3_char 2

//#define L 0
//#define K2_2 2
//#define K2_3 4
//#define K2_3_char 1

//#define L 1
//#define K2_2 4
//#define K2_3 16
//#define K2_3_char 2
//
#define L 2
#define K2_2 8
#define K2_3 64
#define K2_3_char 8

//#define L 3
//#define K2_2 16
//#define K2_3 256
//#define K2_3_char 32

typedef struct node
{
    char * data;
    struct node** child;
    //struct node* next;
    //struct node* parent;    /*just to show that tree is being created correctly*/
}NODE;



typedef struct matrixRep
{
    bitRankW32Int * bt;		//Bitmap representando el arbol
    bitRankW32Int * bn; 	//Bitmap representando las hojas con informacion
    uint bt_len;				//Numero de bits de BT
    uint bn_len;				//Numero de bits de BN
    uint nleaves;				//Numero de hojas con informacion
//    uint K1;
//    uint K2;
//    uint maxRealLevel1; //for divlevel
//    uint maxLevel1;			//Nivel maximo del arbol ref K1
//    uint maxLevel2;			//Nivel maximo del arbol ref K2
    uint cutBt;
    uint lastBt1_len;
    uint numberOfNodes;
    ulong numberOfEdges;
//    uint * div_level_table1;
//    uint * div_level_table2;
    uint * leavesInf;		//Informacion de las hojas
       	FTRep* compressIL;

//    uint * info;
//    uint * element;
//    uint * basex;
//    uint * basey;
//    int iniq;
//    int finq;
    //struct node* next;
    //struct node* parent;    /*just to show that tree is being created correctly*/
}MREP;



typedef struct treeRep{
	  uint part;
    uint tamSubm;
    uint numberOfNodes;
    ulong numberOfEdges;
 //   uint maxLevel;			//Nivel maximo del arbol ref K2
    uint repK1;
    uint repK2;
    
    uint maxRealLevel1; //for divlevel
    uint maxLevel1;			//Nivel maximo del arbol ref K1
    uint maxLevel2;			//Nivel maximo del arbol ref K2

    
    uint * div_level_table1;
    uint * div_level_table2;

    uint * info;
    uint * info2[2];
		
		uint fila;
		uint columna;
		
    uint * element;
    uint * basex;
    uint * basey;

    uint * basep1;
    uint * basep2;
    uint * baseq1;
    uint * baseq2;

    int iniq;
    int finq;
    
    MREP *** submatrices;
 
   	unsigned char * words; //Palabras del vocabulario de hojas ordenadas por frecuencia
  	uint lenWords; //tamaño de las palabras del vocabulario (en bytes)
  	uint zeroNode; //Numero de palabras diferentes del vocabulario de hojas
 
	
} TREP;



#define MAX 100
#ifndef uint
#define uint unsigned int
#endif


uint numberNodes;
uint numberLeaves;
uint numberTotalLeaves;
int max_real_level1, max_Level1, max_Level2;


NODE * createKTree(int _K1, int _K2, int max_real_level1, int maxlevel1, int maxlevel2);

void insertNode(NODE * root, int x, int y);

MREP * createRepresentation(NODE *root, uint numberOfNodes, ulong numberOfEdges);

MREP * loadRepresentation(char * basename);

uint * compactAdjacencyList(TREP * trep, MREP * rep, int x);


uint * compactInverseList(TREP *trep,MREP * rep, int y);




void destroyRepresentation(MREP * rep);

TREP * createTreeRep(uint nodesOrig,ulong edges,uint part,uint subm, uint max_real_level1, uint max_level1, uint max_level2, uint _K1, uint _K2);

void insertIntoTreeRep(TREP * trep, MREP * rep, uint i, uint j);



uint * compactTreeAdjacencyList(TREP * trep, int x);
uint * compactTreeInverseList(TREP * trep, int x);

uint ** compactTreeRangeQuery(TREP * trep, uint p1, uint p2, uint q1, uint q2);
uint compactTreeCheckLink(TREP * trep, uint x, uint y);


void destroyTreeRepresentation(TREP * trep);

void saveTreeRep(TREP * trep, char * basename);
TREP * loadTreeRepresentation(char * basename);

void   compressInformationLeaves(TREP * trep);

void   saveBeforeCompressInformationLeaves(TREP * trep, char * basename);

int comparaFrecListaDesc(const void *a, const void *b);



void  openPartialFile(TREP * trep, char * basename);

void  closePartialFile();

void   partialSave(TREP * trep, MREP * rep, uint fila, uint columna);
void partialdestroyTreeRepresentation(TREP * trep);

#endif
