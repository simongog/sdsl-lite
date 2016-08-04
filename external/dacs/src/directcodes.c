/* Wavelet Tree over Dense Code. -- 
   A word index using wavelet tree strategy over compressed text.

   Programmed by Susana Ladra.

   Author's contact: Susana Ladra, Databases Lab, University of
   A Coruna. Campus de Elvina s/n. Spain  sladra@udc.es

   Copyright (C) 2007  Nieves R. Brisaboa, Antonio Farina and Susana Ladra
   Author's contact: susanaladra@gmail.com

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
  Hash.c: Definition of a HashTable that uses linear hash
  ------------------------------------------------------------------------*/

#include "../include/directcodes.h"

//#define BASE 16
//#define BASE_BITS 4
//#define LONGTABLEBASE 7

//#define REPBSIZE 128
//#define rep->tamtablebase 7
//const unsigned short TABLESUM4[256] = {
//	0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
//1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
//2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
//3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
//4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
//5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
//6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
//7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
//8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
//9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
//10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
//11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
//12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
//13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
//14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
//15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
//
//const byte TABLESUM2[256]={
//	0, 1, 2, 3, 1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6,
//1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6, 4, 5, 6, 7,
//2, 3, 4, 5, 3, 4, 5, 6, 4, 5, 6, 7, 5, 6, 7, 8,
//3, 4, 5, 6, 4, 5, 6, 7, 5, 6, 7, 8, 6, 7, 8, 9,
//1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6, 4, 5, 6, 7,
//2, 3, 4, 5, 3, 4, 5, 6, 4, 5, 6, 7, 5, 6, 7, 8,
//3, 4, 5, 6, 4, 5, 6, 7, 5, 6, 7, 8, 6, 7, 8, 9,
//4, 5, 6, 7, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 10,
//2, 3, 4, 5, 3, 4, 5, 6, 4, 5, 6, 7, 5, 6, 7, 8,
//3, 4, 5, 6, 4, 5, 6, 7, 5, 6, 7, 8, 6, 7, 8, 9,
//4, 5, 6, 7, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 10,
//5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 10, 8, 9, 10, 11,
//3, 4, 5, 6, 4, 5, 6, 7, 5, 6, 7, 8, 6, 7, 8, 9,
//4, 5, 6, 7, 5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 10,
//5, 6, 7, 8, 6, 7, 8, 9, 7, 8, 9, 10, 8, 9, 10, 11,
//6, 7, 8, 9, 7, 8, 9, 10, 8, 9, 10, 11, 9, 10, 11, 12};


//0-8
//const unsigned int rep->tablebase[rep->tamtablebase] = {0,1,257,65793,16843009,4311810305,1103823438081};
//0-2-2-4-8
//const unsigned int rep->tablebase[rep->tamtablebase] = {0,1,5,21,277,65813,16843029};
//0-2-4-8
//const unsigned int rep->tablebase[rep->tamtablebase] = {0,1,5,69,16453,4210757,1077952581};
//0-4-8
//const unsigned int rep->tablebase[rep->tamtablebase] = {0,1,17,4113,1052689,269488145};
//const unsigned int rep->tablebase[7] = {0,16,272,4368,69904,1118480,17895696};
//const unsigned int pot256[4] = {1,256,65536,16777216};

/*-----------------------------------------------------------------

  ---------------------------------------------------------------- */

#define FACT_RANK 20

uint mypow(uint base, uint exponente) {
  int i=0;
  uint result=1;
  for (i=0; i < exponente; i++){
    result *= base;
  }
  return result; 
}




ushort * optimizationk(uint * acumFreqs,int maxInt, uint * nkvalues) {
  int sizeVoc = maxInt;

  //uint listLength = acumFreqs[sizeVoc];	
  uint nBits = bits(sizeVoc);

  //fprintf(stderr,"Length of the list: %d, max bits: %d\n",listLength,nBits);

  //Esta tabla tiene el tamanho que ocupa la mejor opcion para los primeros x bits
  long * tableSize = (long *) malloc(sizeof(long)*(nBits+1));

  //Estos dos arrays indican los bs optimos hasta este punto
  //El primer array dice cuantos niveles hay hasta ese punto y el otro los bits en cada nivel
  uint * tableNLevels = (uint *) malloc(sizeof(uint)*(nBits+1));
  uint ** tableKvalues = (uint **) malloc(sizeof(uint*)*(nBits+1));

  int j;
  ulong maxSize = 0, maxPos = 0;
  int posVocInf, posVocSup;
  ulong currentSize;
  //Para la optimizacion, hay que mirar el máximo de todas las opciones anteriores anhadiendo
  //hasta el bit actual

  tableSize[0]=0;
  tableNLevels[0]=0;
  tableKvalues[0]=NULL;


  int i;
  for(i=1;i<=nBits;i++){
    maxSize = -1;
    maxPos=0;


    for(j=0;j<i;j++){
      if(i==nBits)
        posVocInf=0;
      else
        posVocInf=1<<(nBits-i);
      posVocSup=(1<<(nBits-j));
      if(posVocSup>=sizeVoc)
        posVocSup=sizeVoc;
      if(j==0)
        currentSize = tableSize[j]+	((ulong)(acumFreqs[sizeVoc]-acumFreqs[posVocInf]))*((i-j));	
      else
        currentSize = tableSize[j]+	((ulong)(acumFreqs[sizeVoc]-acumFreqs[posVocInf]))*((i-j)+1)+(acumFreqs[sizeVoc]-acumFreqs[posVocInf])/FACT_RANK;	

      if(maxSize>currentSize){
        maxSize=currentSize;
        maxPos = j;
      }

    }

    tableSize[i] = maxSize;
    tableNLevels[i]=tableNLevels[maxPos]+1;
    tableKvalues[i]=(uint *)malloc(sizeof(uint)*tableNLevels[i]);
    for(j=0; j < tableNLevels[i]-1; j++)
      tableKvalues[i][j]=tableKvalues[maxPos][j];
    tableKvalues[i][tableNLevels[i]-1] = i-maxPos;			


  }


  ulong sumaTotal=0;
  int bitCountInf=0, bitCountSup=0, bitsCount;
  for(j=0;j<tableNLevels[nBits];j++){
    bitsCount = tableKvalues[nBits][tableNLevels[nBits]-1-j];
    bitCountSup+=bitsCount;
    if(bitCountInf==0)
      posVocInf=0;
    else
      posVocInf=1<<bitCountInf;
    posVocSup=(1<<bitCountSup);
    if(posVocSup>=sizeVoc)
      posVocSup=sizeVoc;
    if(j==tableNLevels[nBits])
      sumaTotal += 	((ulong)(acumFreqs[sizeVoc]-acumFreqs[posVocInf]))*bitsCount;	

    else
      sumaTotal += ((ulong)(acumFreqs[sizeVoc]-acumFreqs[posVocInf]))*(bitsCount+1)+(acumFreqs[sizeVoc]-acumFreqs[posVocInf])/FACT_RANK;	




    bitCountInf+=bitsCount;
  }



  (*nkvalues)=tableNLevels[nBits];

  ushort * kvalues = (ushort *)malloc(sizeof(ushort)*tableNLevels[nBits]);

  for(j=0;j<tableNLevels[nBits];j++){
    bitsCount = tableKvalues[nBits][tableNLevels[nBits]-1-j];

    kvalues[j] = bitsCount;



    bitCountInf+=bitsCount;
  }


  free (tableSize);
  for (i = 1; i <= nBits; ++i)
    free (tableKvalues[i]);
  free (tableKvalues);
  free (tableNLevels);
  return kvalues;
}

FTRep* createFT(uint *list,uint listLength){
  FTRep * rep = (FTRep *) malloc(sizeof(struct sFTRep));
  uint *levelSizeAux;
  uint *cont;	
  uint *contB;

  ushort* kvalues;
  uint nkvalues;

  rep->listLength = listLength;
  register uint i;
  int j, k;
  uint l;
  uint value, newvalue;
  uint bits_BS_len = 0;

  //ushort kvalues[4] = {0,2,4,8};
  //uint nkvalues=4;

  uint maxInt=0;

  for(i=0; i < listLength; i++)
    if(maxInt < list[i])
      maxInt = list[i];

  //fprintf(stderr,"maxInt : %d\n",maxInt);

  maxInt++;

  uint * weight = (uint *) malloc(sizeof(uint)*maxInt);


  for(l=0; l<maxInt; l++)
    weight[l] = 0;

  for(i=0;i<listLength;i++)
    weight[list[i]]++;


  uint * acumFreq = (uint *) malloc(sizeof(uint)*(maxInt+1));

  acumFreq[0]=0;	
  for(i=0;i<maxInt;i++)
    acumFreq[i+1]=acumFreq[i]+weight[i];

  free(weight);


  kvalues = optimizationk(acumFreq,maxInt,&nkvalues);


  free(acumFreq);



  ushort kval;
  uint oldval =0;
  uint newval =0;

  i=0;
  uint multval=1;
  do{
    oldval=newval;
    if(i>=nkvalues){
      kval = 1<<(kvalues[nkvalues-1]);
    }
    else
      kval=1<<(kvalues[i]);	
    multval*=kval;
    newval = oldval+multval;//mypow(kval,i);
    i++;
  }
  while(oldval<newval);

  rep->tamtablebase = i;
  rep->tablebase = (uint *) malloc(sizeof(uint)*rep->tamtablebase);
  levelSizeAux = (uint *) malloc(sizeof(uint)*rep->tamtablebase);
  cont = (uint *) malloc(sizeof(uint)*rep->tamtablebase);
  contB = (uint *) malloc(sizeof(uint)*rep->tamtablebase);

  oldval =0;
  newval =0;
  multval=1;	
  for(i=0;i<rep->tamtablebase;i++){
    oldval=newval;
    if(i>=nkvalues){
      kval = 1<<(kvalues[nkvalues-1]);
    }
    else
      kval=1<<(kvalues[i]);	
    multval*=kval;
    newval = oldval+multval;//mypow(kval,i);
    rep->tablebase[i]=oldval;
  }	


  for(i=0;i<rep->tamtablebase;i++){
    levelSizeAux[i]=0;

  }
  //Reservando espacio para los niveles

  for (i=0;i<listLength;i++){
    value = list[i];
    for(j=0;j<rep->tamtablebase;j++)
      if(value>=rep->tablebase[j])
        levelSizeAux[j]++;
  }
  //
  j=0;
  //	//Contadores
  while((j<rep->tamtablebase)&&(levelSizeAux[j]!=0)){
    j++;
  }
  rep->nLevels = j;
  rep->levelsIndex = (uint *) malloc(sizeof(uint)*(rep->nLevels+1));
  bits_BS_len =0;

  rep->base = (uint *)malloc(sizeof(uint)*rep->nLevels);
  rep->base_bits = (ushort *)malloc(sizeof(ushort)*rep->nLevels);

  for(i=0;i<rep->nLevels;i++){
    if(i>=nkvalues){
      rep->base[i]=1<<(kvalues[nkvalues-1]);
      rep->base_bits[i]=kvalues[nkvalues-1];
    }
    else{
      rep->base[i]=1<<(kvalues[i]);
      rep->base_bits[i]=kvalues[i];
    }
  }

  uint tamLevels =0;



  tamLevels=0;
  for(i=0;i<rep->nLevels;i++)
    tamLevels+=rep->base_bits[i]*levelSizeAux[i];

  rep->iniLevel = (uint *)malloc(sizeof(uint)*rep->nLevels);		
  rep->tamCode=tamLevels;
  uint indexLevel=0;
  rep->levelsIndex[0]=0;
  for(j=0;j<rep->nLevels;j++){
    rep->levelsIndex[j+1]=rep->levelsIndex[j] + levelSizeAux[j];
    //if(j>0){
    rep->iniLevel[j] = indexLevel;
    cont[j]=rep->iniLevel[j];
    indexLevel+=levelSizeAux[j]*rep->base_bits[j];
    //}
    contB[j]=rep->levelsIndex[j];

  }


  rep->levels = (uint *) malloc(sizeof(uint)*(tamLevels/W+1));
  //fprintf(stderr,"tamaño de levels: %d\n",sizeof(uint)*(tamLevels/W+1));
  bits_BS_len = rep->levelsIndex[rep->nLevels-1]+1; 
  //Se pone el ultimo a 0 para ahorrarnos comparaciones despues en la descompresion
  uint * bits_BS = (uint *) malloc(sizeof(uint)*(bits_BS_len/W+1));
  for(i=0; i<((bits_BS_len)/W+1);i++)
    bits_BS[i]=0;
  for(i=0;i<listLength;i++){
    value = list[i];
    j=rep->nLevels-1;

    while(j>=0){
      if(value >= rep->tablebase[j]){

        newvalue = value- rep->tablebase[j];

        for(k=0;k<j;k++){



          bitwrite(rep->levels,cont[k],rep->base_bits[k],newvalue%rep->base[k]);
          cont[k]+=rep->base_bits[k];
          contB[k]++;
          newvalue = newvalue/rep->base[k];
        }
        k=j;


        bitwrite(rep->levels,cont[j],rep->base_bits[j],newvalue%rep->base[j]);
        cont[j]+=rep->base_bits[j];
        contB[j]++;
        if(j<rep->nLevels-1){
          bitset(bits_BS,contB[j]-1);
        }

        break;
      }
      j--;
    }
    //Para j=0 solo se cubre el bitmap

  }
  //Para simular ultimo array:

  bitset(bits_BS,bits_BS_len-1);
  //rep->bits_bitmap = bits_BS;
  rep->bS = createBitRankW32Int(bits_BS, bits_BS_len , 1, 20); 	




  rep->rankLevels = (uint *) malloc(sizeof(uint)*rep->nLevels);
  for(j=0;j<rep->nLevels;j++)
    rep->rankLevels[j]= rank(rep->bS, rep->levelsIndex[j]-1);


  free(cont);
  free(contB);
  free(levelSizeAux);
  free(kvalues);
  return rep;
}


/*-----------------------------------------------------------------

  ---------------------------------------------------------------- */

uint accessFT(FTRep * listRep,uint param){
  uint mult=0;
  //register uint i;
  register uint j;
  //uint * rankLevel = listRep->rankLevels;
  //byte * list;
  //uint n , partialSum=0, sumAux=0;
  uint partialSum=0;
  //uint ini = param-1;
  uint ini = param;
  //bitRankW32Int * bS = listRep->bS;
  //uint * bsData = listRep->bS->data;
  uint nLevels=listRep->nLevels;
  //uint levelIndex;
  uint * level;
  uint readByte;
  uint cont,pos, rankini;


  //	fprintf(stderr,"Queriendo leer la posicion: %d\n",ini);
  partialSum=0;
  j=0;
  level=listRep->levels ;//+ (listRep->levelsIndex[j]>>1);
  //		fprintf(stderr,"primera posicion del array levels: %d %d\n",level[0],(byte)level[0]);
  //cont=ini+(listRep->levelsIndex[j]&0x1);
  pos=listRep->levelsIndex[j]+ini;

  mult=0;
  //readByte = ((*(level+(cont>>1)))>>(BASE_BITS*(cont&0x1)))&0xF;
  cont = listRep->iniLevel[j]+ini*listRep->base_bits[j];
  //		fprintf(stderr,"leyendo los %d bits \n",listRep->base_bits[j]);
  //		fprintf(stderr,"inilevel: %d \n",listRep->iniLevel[j]);
  //		fprintf(stderr,"ini: %d\n",ini);
  //		fprintf(stderr,"base_bits: %d\n",listRep->base_bits[j]);		
  //
  //		fprintf(stderr,"%d posicion\n",cont);
  //		fprintf(stderr,"nivel %d\n",j);		
  //		fprintf(stderr,"leyendo los %d bits desde %d posicion del nivel %d\n",listRep->base_bits[j],cont,j);
  readByte = bitread(level,cont,listRep->base_bits[j]);
  //		fprintf(stderr,"readbyte: %d\n",readByte);
  //fprintf(stderr,"contenido del vector: %d readByte: %d con cont[%d]: %d\n",*(level[j]+(cont[j]>>1)),readByte,j,cont[j]);
  //fprintf(stderr,"readByte: %d... %d\n",readByte,(*(level[j]+(cont[j]>>1))));
  //	fprintf(stderr,"readByte: %d\n",readByte);
  //fprintf(stderr,"pos[%d]= %d\n",j,pos[j]);
  if(nLevels == 1){
    return readByte;
  }
  while((!bitget(listRep->bS->data,pos))){
    //fprintf(stderr,"pos[%d]= %d\n",j,pos);
    rankini = rank(listRep->bS, listRep->levelsIndex[j]+ini-1) - listRep->rankLevels[j];
    ini = ini-rankini;

    //fprintf(stderr,"readByte: %d\n",readByte);
    //partialSum = partialSum+ readByte*pot256[j];
    partialSum = partialSum+ (readByte<<mult);

    mult+=listRep->base_bits[j];
    j++;
    //fprintf(stderr,"pos[%d] = %d lsit-> levels %d cont[%d] = %d \n",j-1,pos[j-1]-1,listRep->levels, j,cont[j]);

    //level=listRep->levels + (listRep->levelsIndex[j]>>1);
    //cont=ini+(listRep->levelsIndex[j]&0x1);
    cont = listRep->iniLevel[j]+ini*listRep->base_bits[j];
    pos=listRep->levelsIndex[j]+ini;

    //readByte = ((*(level+(cont>>1)))>>(BASE_BITS*(cont&0x1)))&0xF;
    //		fprintf(stderr,"leyendo los %d bits desde %d posicion del nivel %d\n",listRep->base_bits[j],cont,j);
    readByte = bitread(level,cont,listRep->base_bits[j]);
    //		fprintf(stderr,"readbyte: %d\n",readByte);
    //fprintf(stderr,"contenido del vector: %d readByte: %d con cont[%d]: %d\n",*(level[j]+(cont[j]>>1)),readByte,j,cont[j]);

    //fprintf(stderr,"readByte: %d... %d\n",readByte,(*(level+(cont>>1))));
    if(j==nLevels-1){
      break;
    }


  }
  //fprintf(stderr,"readByte: %d\n",readByte);
  //partialSum = partialSum + readByte*pot256[j] + rep->tablebase[j];
  partialSum = partialSum + (readByte<<mult) + listRep->tablebase[j];
  //fprintf(stderr,"partialSum = %d, param = %u \n",partialSum,param);
  //list[i]=sum
  //	fprintf(stderr,"devolviendo %d\n",partialSum);
  return partialSum;

}

bool equalsRank(bitRankW32Int *lhs, bitRankW32Int *rhs) {
  if (lhs->owner != rhs->owner || lhs->integers != rhs->integers ||
      lhs->factor != rhs->factor || lhs->b != rhs->b ||
      lhs->s != rhs->s || lhs->n != rhs->n)
    return false;

  uint i;
  for (i = 0; i < lhs->n/W + 1; ++i)
    if (lhs->data[i] != rhs->data[i])
      return false;

  for (i = 0; i < lhs->n/lhs->s + 1; ++i)
    if (lhs->data[i] != rhs->data[i])
      return false;

  return true;
}
bool equalsFT(FTRep *lhs, FTRep *rhs) {
  if (lhs->listLength != rhs->listLength || lhs->nLevels != rhs->nLevels ||
      lhs->tamCode != rhs->tamCode || lhs->tamtablebase != rhs->tamtablebase)
    return false;

  uint i;
  for (i = 0; i < lhs->nLevels; ++i) {
    if (lhs->base_bits[i] != rhs->base_bits[i] ||
        lhs->base[i] != rhs->base[i] ||
        lhs->iniLevel[i] != rhs->iniLevel[i] ||
        lhs->rankLevels[i] != rhs->rankLevels[i] ||
        lhs->levelsIndex[i] != rhs->levelsIndex[i])
      return false;
  }

  if (lhs->levelsIndex[lhs->nLevels] != rhs->levelsIndex[lhs->nLevels])
    return false;

  for (i = 0; i< lhs->tamtablebase; ++i)
    if (lhs->tablebase[i] != rhs->tablebase[i])
      return false;

  for (i = 0; i < lhs->tamCode/W + 1; ++i)
    if (lhs->levels[i] != rhs->levels[i])
      return false;
  return equalsRank(lhs->bS, rhs->bS);
}



/*-----------------------------------------------------------------

  ---------------------------------------------------------------- */

uint * decompressFT(FTRep * listRep, uint n){
  uint mult=0;
  register uint i;
  register uint j;
  //uint partialSum=0, sumAux=0;
  uint partialSum=0;
  //uint ini = 0;
  //bitRankW32Int * bS = listRep->bS;
  //uint * bsData = listRep->bS->data;
  uint nLevels=listRep->nLevels;
  //uint levelIndex;
  uint * level=listRep->levels;
  byte readByte;
  uint * list = (uint *) malloc(sizeof(uint)*n);
  uint * cont = (uint *) malloc(sizeof(byte*)*listRep->nLevels);
  uint * pos = (uint *) malloc(sizeof(uint)*listRep->nLevels);

  for(j=0;j<nLevels;j++){
    cont[j]=listRep->iniLevel[j];
    pos[j]=listRep->levelsIndex[j];
    //fprintf(stderr,"cont[%d] = %d, pos[%d]= %d\n",j,cont[j],j,pos[j]);
  }
  // uint * rankLevel = listRep->rankLevels;
  // byte * list;

  //uint cont,pos, rankini;
  //fprintf(stderr,"Queriendo leer la posicion: %d\n",ini);
  for(i=0;i<n;i++){
    partialSum=0;
    j=0;
    //level=listRep->levels + (listRep->levelsIndex[j]>>1);
    //cont=ini+(listRep->levelsIndex[j]&0x1);
    //pos=listRep->levelsIndex[j]+ini;

    mult=0;
    readByte = bitread(level,cont[j],listRep->base_bits[j]);
    cont[j]+=listRep->base_bits[j];
    //readByte = ((*(level+(cont[j]>>1)))>>(BASE_BITS*(cont[j]&0x1)))&0xF;
    //fprintf(stderr,"contenido del vector: %d readByte: %d con cont[%d]: %d\n",*(level[j]+(cont[j]>>1)),readByte,j,cont[j]);
    //fprintf(stderr,"readByte: %d... %d\n",readByte,(*(level[j]+(cont[j]>>1))));

    //fprintf(stderr,"pos[%d]= %d\n",j,pos[j]);
    while((!bitget(listRep->bS->data,pos[j]))){
      //fprintf(stderr,"pos[%d]= %d\n",j,pos[j]);
      //rankini = rank(listRep->bS, listRep->levelsIndex[j]+ini-1) - listRep->rankLevels[j];
      //ini = ini-rankini;
      pos[j]++;
      //cont[j]++;
      //fprintf(stderr,"readByte: %d\n",readByte);
      //partialSum = partialSum+ readByte*pot256[j];
      partialSum = partialSum+ (readByte<<mult);
      mult+=listRep->base_bits[j];
      j++;

      //fprintf(stderr,"pos[%d] = %d lsit-> levels %d cont[%d] = %d \n",j-1,pos[j-1]-1,listRep->levels, j,cont[j]);


      //level=listRep->levels + (listRep->levelsIndex[j]>>1);
      //cont=ini+(listRep->levelsIndex[j]&0x1);
      //pos=listRep->levelsIndex[j]+ini;
      readByte = bitread(level,cont[j],listRep->base_bits[j]);
      cont[j]+=listRep->base_bits[j];
      //readByte = ((*(level+(cont[j]>>1)))>>(BASE_BITS*(cont[j]&0x1)))&0xF;
      //fprintf(stderr,"contenido del vector: %d readByte: %d con cont[%d]: %d\n",*(level[j]+(cont[j]>>1)),readByte,j,cont[j]);

      //fprintf(stderr,"readByte: %d... %d\n",readByte,(*(level[j]+(cont[j]>>1))));
      if(j==nLevels-1){
        break;
      }

    }
    //fprintf(stderr,"readByte: %d\n",readByte);
    //partialSum = partialSum + readByte*pot256[j] + rep->tablebase[j];
    if(j<nLevels-1){
      pos[j]++;
    }


    partialSum = partialSum + (readByte<<mult) + listRep->tablebase[j];
    //fprintf(stderr,"sum = %d, partialSum = %d, param = %u \n",sum,partialSum,param);
    //list[i]=sum
    //if(i==16404)
    //fprintf(stderr,"devolviendo %d\n",partialSum);
    list[i]=partialSum;
  }
  free(cont);
  free(pos);
  return list;
}


/*-----------------------------------------------------------------

  ---------------------------------------------------------------- */
void destroyFT(FTRep * rep){
  free(rep->levelsIndex);
  free(rep->iniLevel);
  free(rep->rankLevels);
  free(rep->tablebase);
  free(rep->levels);
  //destroyBitmap(rep->bS);
  destroyBitRankW32Int(rep->bS);
  //free(rep->bits_bitmap);
  free(rep->base);
  free(rep->base_bits);
  free(rep);
}

/*-----------------------------------------------------------------

  ---------------------------------------------------------------- */
void saveFT(FTRep * rep, FILE * flist){
  //int i;

  fwrite(&(rep->listLength),sizeof(uint),1,flist);
  fwrite(&(rep->nLevels),sizeof(byte),1,flist);
  fwrite(&(rep->tamCode),sizeof(uint),1,flist);
  fwrite(&(rep->tamtablebase),sizeof(uint),1,flist);
  fwrite(rep->tablebase,sizeof(uint),rep->tamtablebase,flist);	
  fwrite(rep->base_bits,sizeof(ushort),rep->nLevels,flist);
  fwrite(rep->base,sizeof(uint),rep->nLevels,flist);
  fwrite(rep->levelsIndex,sizeof(uint),rep->nLevels+1,flist);
  fwrite(rep->iniLevel,sizeof(uint),rep->nLevels,flist);
  fwrite(rep->rankLevels,sizeof(uint),rep->nLevels,flist);
  //for(i=0;i<rep->nLevels;i++)
  //	fprintf(stderr,"i:%d ranklevel: %d\n",i,rep->rankLevels[i]);
  fwrite(rep->levels,sizeof(uint),(rep->tamCode/W+1),flist);
  //fwrite(rep->bS->data,sizeof(uint),((rep->bS->n + W -1)/W) ,flist);	
  //fwrite(rep->bits_bitmap,sizeof(uint),((rep->levelsIndex[rep->nLevels-1]+1 + W -1)/W),flist);
  //	fwrite(rep->bS->data,sizeof(uint),652,flist);
  //fprintf(stderr,"accediendo a data: %x, W:%d, bytes: %d\n",rep->bS->data[rep->bS->n/W],W,rep->bS->n/W+1);
  save(rep->bS,flist);
  //	fprintf(stderr,"accediendo a data: %x, W:%d, bytes: %d\n",rep->bS->data[rep->bS],W,rep->bS->n/W+1);

  //fclose(flist);

  //	char *filename2 = (char *) malloc(sizeof(char)*(strlen(filename)+4));
  //  strcpy(filename2,filename);
  //  strcat(filename2,".bm");
  //  saveBitmap(filename2,rep->bS);


  //	fprintf(stderr,"bitmap: %d, sumas: %d\n",sizeof(uint) * rep->bS->sSize + sizeof(byte) * rep->bS->bSize,sizeof(uint)*((rep->stoppersSize*rep->sumBits + W -1)/W));
  //	uint espSinEstr, espEstrDisco, espEstrMemoria, espEstr, espTotal;
  //	espSinEstr =sizeof(uint)*(rep->nLevels+1) + sizeof(byte)*rep->levelsIndex[rep->nLevels]/2+sizeof(uint) *(rep->bS->n/W+1);
  //	espEstrDisco = sizeof(byte) + sizeof(uint)*7 + sizeof(uint)*((rep->stoppersSize*rep->sumBits + W -1)/W)+sizeof(uint) * rep->bS->sSize + sizeof(byte) * rep->bS->bSize;
  //	espEstrMemoria = sizeof(uint *)*2 + sizeof(byte *) +  sizeof(bitmap)-sizeof(uint)*5;
  //	espEstr = espEstrDisco + espEstrMemoria;
  //	espTotal = espEstr + espSinEstr;
  //fprintf(stderr,"\n---\nSin estructuras ocupa: %d\n",espSinEstr);
  //fprintf(stderr,"Las estructuras adicionales ocupan: %d\n",espEstr);
  //fprintf(stderr,"Espacio total: %d, espacio del fichero: %d\n", espEstrDisco+espSinEstr,fileSize(filename)+fileSize(filename2));
  //fprintf(stderr,"Porcentaje de estructuras: %f sobre el total\n---\n", ((float)espEstr)*100/(4*rep->listLength));


}
FTRep* loadFT(FILE * flist){
  //FILE * flist;

  //int i;
  FTRep * rep = (FTRep *) malloc(sizeof(struct sFTRep));
  //flist = fopen(filename,"r");
  fread(&(rep->listLength),sizeof(uint),1,flist);
  fread(&(rep->nLevels),sizeof(byte),1,flist);
  fread(&(rep->tamCode),sizeof(uint),1,flist);

  fread(&(rep->tamtablebase),sizeof(uint),1,flist);
  rep->tablebase = (uint *) malloc(sizeof(uint)*rep->tamtablebase);
  fread(rep->tablebase,sizeof(uint),rep->tamtablebase,flist);	

  rep->base_bits = (ushort *) malloc(sizeof(ushort)*rep->nLevels);
  fread(rep->base_bits,sizeof(ushort),rep->nLevels,flist);
  //for(i=0;i<rep->nLevels;i++)
  //	fprintf(stderr,"%d,",rep->base_bits[i]);
  //fprintf(stderr," son los bits usados\n");

  rep->base = (uint *) malloc(sizeof(uint)*rep->nLevels);
  fread(rep->base,sizeof(uint),rep->nLevels,flist);
  //for(i=0;i<rep->nLevels;i++)
  //	fprintf(stderr,"base[%d]=%d\n",i,rep->base[i]);

  rep->levelsIndex = (uint *) malloc(sizeof(uint)*(rep->nLevels+1));
  fread(rep->levelsIndex,sizeof(uint),rep->nLevels+1,flist);

  rep->iniLevel = (uint *) malloc(sizeof(uint)*rep->nLevels);
  fread(rep->iniLevel,sizeof(uint),rep->nLevels,flist);

  rep->rankLevels = (uint *) malloc(sizeof(uint)*(rep->nLevels));
  fread(rep->rankLevels,sizeof(uint),rep->nLevels,flist);
  //for(i=0;i<rep->nLevels;i++)
  //	fprintf(stderr,"i:%d ranklevel: %d\n",i,rep->rankLevels[i]);
  rep->levels = (uint *) malloc(sizeof(uint)*(rep->tamCode/W+1));	
  fread(rep->levels,sizeof(uint),(rep->tamCode/W+1),flist);

  //rep->bits_bitmap = (uint *) malloc(sizeof(uint)*((rep->levelsIndex[rep->nLevels-1]+1 + W -1)/W));
  //fread(rep->bits_bitmap,sizeof(uint),((rep->levelsIndex[rep->nLevels-1]+1 + W -1)/W),flist);

  rep->bS = (bitRankW32Int *) malloc(sizeof(struct sbitRankW32Int));
  load(rep->bS,flist);	

  //fclose(flist);

  //fprintf(stderr,"levels:%d\n",rep->nLevels);
  //	fprintf(stderr,"blockSize: %d. stoppersSize: %d, sumBits:%d\n",rep->blockSize,rep->stoppersSize,rep->sumBits);
  //	char *filename2 = (char *) malloc(sizeof(char)*(strlen(filename)+4));
  //  strcpy(filename2,filename);
  //  strcat(filename2,".bm");
  //  rep->bS = loadBitmap(filename2,rep->bits_bitmap,rep->levelsIndex[rep->nLevels-1]+1);



  //	uint espSinEstr, espEstrDisco, espEstrMemoria, espEstr, espTotal;
  //	espSinEstr =sizeof(uint)*(rep->nLevels+1) + sizeof(byte)*rep->levelsIndex[rep->nLevels]/2+sizeof(uint) *(rep->bS->n/W+1);
  //	espEstrDisco = sizeof(byte) + sizeof(uint)*7 + sizeof(uint)*((rep->stoppersSize*rep->sumBits + W -1)/W)+sizeof(uint) * rep->bS->sSize + sizeof(byte) * rep->bS->bSize;
  //	espEstrMemoria = sizeof(uint *)*2 + sizeof(byte *) +  sizeof(bitmap)-sizeof(uint)*5;
  //	espEstr = espEstrDisco + espEstrMemoria;
  //	espTotal = espEstr + espSinEstr;
  //fprintf(stderr,"\n---\nSin estructuras ocupa: %d\n",espSinEstr);
  //fprintf(stderr,"Las estructuras adicionales ocupan: %d\n",espEstr);
  //fprintf(stderr,"Espacio total: %d, espacio del fichero: %d\n", espEstrDisco+espSinEstr,fileSize(filename)+fileSize(filename2));
  //fprintf(stderr,"Porcentaje de estructuras: %f sobre el total\n---\n", ((float)espEstr)*100/(4*rep->listLength));


  return rep;
}


uint memoryUsage(FTRep* rep) {
  return sizeof(uint)*rep->tamtablebase 
    + sizeof(ushort)*rep->nLevels
    + sizeof(ushort)*rep->nLevels
    + sizeof(uint)*(rep->nLevels+1)
    + sizeof(uint)*(rep->nLevels)
    + sizeof(uint)*(rep->nLevels)
    + sizeof(uint)*(rep->tamCode/W+1)
    + spaceRequirementInBits(rep->bS)/8
    + sizeof(struct sFTRep);
}
