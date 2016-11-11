#include <stdio.h>
#include <math.h>
#include <string.h>
#include "kTree.h"
#include <sys/time.h>
#include <sys/resource.h>

#define swap( x, y ) { int temp; temp=*(x); *(x)=*(y); *(y)=temp; }

/* Time meassuring */
double ticks;
struct tms t1,t2;

void start_clock() {
	times (&t1);
}

double stop_clock() {
	times (&t2);
	return (t2.tms_utime-t1.tms_utime)/ticks;
}

int main(int argc, char* argv[]){
	ticks= (double)sysconf(_SC_CLK_TCK);
  	start_clock();
	FILE *f;
	uint nodes; 
	ulong edges;
	register ulong i;
	
	if(argc<3){
		fprintf(stderr,"USAGE: %s <base name> <hash size>\n",argv[0]);
		return(-1);
	}

	
	//GRABAR IL sin comprimir
	
	uint _part, _tamSubm, _numberOfNodes,_repK1, _repK2,_maxRealLevel1, _maxLevel1,_maxLevel2;
	ulong _numberOfEdges;
	unsigned long addrInTH;
	unsigned int zeroNode;
	unsigned int * positionInTH;

	ulong _totalLeaves;
	ulong hashsize = atol(argv[2]);

	char *filename = (char *) malloc(sizeof(char)*(strlen(argv[1])+5));
	
	strcpy(filename,argv[1]);
	strcat(filename,".il");
	FILE * fvil = fopen(filename,"r");

	fread(&(_part),sizeof(uint),1,fvil);
	fread(&(_tamSubm),sizeof(uint),1,fvil);
	fread(&(_numberOfNodes),sizeof(uint),1,fvil);
	fread(&(_numberOfEdges),sizeof(ulong),1,fvil);	  
	fread(&(_repK1),sizeof(uint),1,fvil);
	fread(&(_repK2),sizeof(uint),1,fvil);
	fread(&(_maxRealLevel1),sizeof(uint),1,fvil);
	fread(&(_maxLevel1),sizeof(uint),1,fvil);
	fread(&(_maxLevel2),sizeof(uint),1,fvil);
	uint **_partnumberOfNodes, **_partcutBt, **_partnleaves, **_partlastBt1_len, ***_partleavesInf;
	ulong ** _partnumberOfEdges;
	FTRep*** _partcompressIL;



	/* SEGUNDA PARTE DE SAVETREE*/
	
	strcpy(filename,argv[1]);
	strcat(filename,".voc");
	FILE * fv = fopen(filename,"w");

	fwrite(&(_part),sizeof(uint),1,fv);
	fwrite(&(_tamSubm),sizeof(uint),1,fv);

	fwrite(&(_numberOfNodes),sizeof(uint),1,fv);
	fwrite(&(_numberOfEdges),sizeof(ulong),1,fv);

	fwrite(&(_repK1),sizeof(uint),1,fv);
	fwrite(&(_repK2),sizeof(uint),1,fv);
	fwrite(&(_maxRealLevel1),sizeof(uint),1,fv);
	fwrite(&(_maxLevel1),sizeof(uint),1,fv);
	fwrite(&(_maxLevel2),sizeof(uint),1,fv);





	_partnumberOfNodes = (uint **)malloc(sizeof(uint*)*_part);
	_partcutBt = (uint **)malloc(sizeof(uint*)*_part);
	_partnleaves = (uint **)malloc(sizeof(uint*)*_part);
	_partlastBt1_len = (uint **)malloc(sizeof(uint*)*_part);
	_partleavesInf = (uint ***)malloc(sizeof(uint**)*_part);
	_partnumberOfEdges = (ulong **)malloc(sizeof(ulong*)*_part);
	_partcompressIL = (FTRep***) malloc(sizeof(FTRep**)*_part);

	for(i=0;i<_part;i++){
		_partnumberOfNodes[i] = (uint *)malloc(sizeof(uint)*_part);
		_partcutBt[i] = (uint *)malloc(sizeof(uint)*_part);
		_partnleaves[i] = (uint *)malloc(sizeof(uint)*_part);
		_partlastBt1_len[i] = (uint *)malloc(sizeof(uint)*_part);
		_partleavesInf[i] = (uint **)malloc(sizeof(uint*)*_part);
		_partnumberOfEdges[i] = (ulong *)malloc(sizeof(ulong)*_part);	
		_partcompressIL[i] = (FTRep**) malloc(sizeof(FTRep*)*_part);		
	}

	uint fila,columna;
	for(fila=0;fila<_part;fila++){
		for(columna=0;columna<_part;columna++){  
			fread(&(_partnumberOfNodes[fila][columna]),sizeof(uint),1,fvil);
			fread(&(_partnumberOfEdges[fila][columna]),sizeof(ulong),1,fvil);

			if(_partnumberOfEdges[fila][columna]==0)
				continue;
			fread(&(_partcutBt[fila][columna]),sizeof(uint),1,fvil);
			fread(&(_partlastBt1_len[fila][columna]),sizeof(uint),1,fvil);
			  //fwrite(&(rep->maxLevel),sizeof(uint),1,fi);
			fread(&(_partnleaves[fila][columna]),sizeof(uint),1,fvil);
			  //fwrite (rep->leavesInf,sizeof(uint),_partnleaves[fila][columna]*K2*K2/W+1,fi);
			_partleavesInf[fila][columna] = (uint *)malloc(sizeof(uint)*_partnleaves[fila][columna]*(K2_2*K2_2+W-1)/W);
			fread(_partleavesInf[fila][columna],sizeof(uint),(_partnleaves[fila][columna]*K2_2*K2_2+W-1)/W,fvil);

		}
	}

	fread(&(_totalLeaves),sizeof(ulong),1,fvil);
	
	
	

	unsigned char * ilchar = (unsigned char *) malloc(sizeof(unsigned char)*_totalLeaves*K2_3_char);
	
	unsigned char * ilchar2 = ilchar;
	uint j;
	
	
	
	
	
	unsigned char *aWord;
	unsigned int size;
	
	ulong contadorTL = 0;
	for(fila=0;fila<_part;fila++){
		for(columna=0;columna<_part;columna++){
			if(_partnumberOfEdges[fila][columna]==0)
				continue;
			for(i=0;i<_partnleaves[fila][columna];i++){
				//Usando una hash para el vocabulario...
				aWord=ilchar2;  //the word parsed.
				
				for(j=0;j<K2_3;j++){
					if(bitget(_partleavesInf[fila][columna],i*K2_3+j))
						bitsetchar(aWord,j);
					else
						bitcleanchar(aWord,j);
					
				}
				size= K2_3_char;
				//Guardando
				ilchar2+=K2_3_char;
				contadorTL++;
				
			}
		}
	}

	for(fila=0;fila<_part;fila++){
		for(columna=0;columna<_part;columna++){

			free(_partleavesInf[fila][columna]);
		}

	}
	uint lenWords = K2_3_char;
	
	zeroNode = 0;
	
	_memMgr = createMemoryManager();

	initialize_hash(hashsize);

	positionInTH = (unsigned int*) malloc (_totalLeaves * sizeof(unsigned int));

	unsigned char * lastpost= &(ilchar[_totalLeaves*K2_3_char]); 
	//Creación del vocabulario
	ulong ilpos = 0;
	for(i=0;i<_totalLeaves;i++){
		//Usando una hash para el vocabulario...
		aWord=&(ilchar[ilpos]);  //the word parsed.

		size= K2_3_char;
		j = search ((unsigned char *)aWord, size, &addrInTH );

		if (j==zeroNode) {
			insertElement ((unsigned char *) aWord, size, &addrInTH);
			hash[addrInTH].weight = 0;
			hash[addrInTH].size = 0;
			hash[addrInTH].len = K2_3_char;
					//printf("%d de %d\n",zeroNode,_totalLeaves);
			positionInTH[zeroNode] = addrInTH;
			zeroNode++;
		}

		hash[addrInTH].weight +=1;

		ilpos+=K2_3_char;		


	}		

	//printf("Zeronode: %d, ilpos: %d\n",zeroNode, ilpos);



	//Compresion de hojas

	int k=0;
	// Sorting the vocabulary by frequency.

	{	//Moves all the words with frequency = 1 to the end of the vocabulary.
		register int ii;
		register int kk;

		kk=zeroNode-1;
		ii=0;
		while (ii<kk){
			while ((hash[positionInTH[ii]].weight!=1) && (ii <kk)) { ii++; }
			while ((hash[positionInTH[kk]].weight==1) && (ii <kk)) { kk--; }

			if (ii<kk){
				swap(&positionInTH[ii], &positionInTH[kk]);
				kk--;
				ii++;
			}
		}

		//k=ii; 
		k=ii+1; //the lenght of the vector to be sorted with qsort. So v[0 .. k-1]
	}

	//printf("Aplies qsort to the words with frequency > 1.\n");
	//Aplies qsort to the words with frequency > 1.
	qsort(positionInTH,k,sizeof(unsigned int),comparaFrecListaDesc);

	int jj;
	//Generates codes sequentially

	
	uint totalLeavesCount=0;
	for(i=0;i<zeroNode;i++){
		hash[positionInTH[i]].codeword = i;
		totalLeavesCount += hash[positionInTH[i]].weight;
	}
	//fprintf(stderr,"totalLeavesCount: %d, numberLeaves: %d\n",totalLeavesCount, numberLeaves);
	

	/********************** Beginning of the second pass **********************/
	
	uint icont=0;

	unsigned int codeword;
	unsigned int tam;
	

	//Compactando información de las hojas...

	ilpos=0;
	for(fila=0;fila<_part;fila++){
		for(columna=0;columna<_part;columna++){
			if(_partnumberOfEdges[fila][columna]==0)
				continue;
			unsigned int tamTotal = 0;
			uint * listIL = (uint *) malloc(sizeof(uint)*_partnleaves[fila][columna]);
			uint listILCount =0;


			for(i=0;i<_partnleaves[fila][columna];i++){


				aWord=&(ilchar[ilpos]);  //the word parsed.
				size = K2_3_char;
				j = search ((unsigned char *)aWord, size, &addrInTH );
						
				listIL[listILCount++]=hash[addrInTH].codeword;

				ilpos+=K2_3_char;		
			}



			_partcompressIL[fila][columna] = createFT(listIL,_partnleaves[fila][columna]);

			free(listIL);
		}
	}

   	unsigned char * words; //Palabras del vocabulario de hojas ordenadas por frecuencia

   	words = (unsigned char *) malloc(sizeof(unsigned char)*zeroNode*lenWords);

   	int wc = 0;
   	for (i=0;i<zeroNode;i++){
   		for(j=0;j<lenWords;j++){
   			words[wc++]=hash[positionInTH[i]].word[j];
   		}
   	}

   	free(ilchar);


	/* SEGUNDA PARTE DE SAVETREE*/
/*	
			strcpy(filename,argv[1]);
		  strcat(filename,".voc");
		  FILE * fv = fopen(filename,"w");
		
		  fwrite(&(_part),sizeof(uint),1,fv);
		  fwrite(&(_tamSubm),sizeof(uint),1,fv);
		
		  fwrite(&(_numberOfNodes),sizeof(uint),1,fv);
		  fwrite(&(_numberOfEdges),sizeof(ulong),1,fv);
		  
		  fwrite(&(_repK1),sizeof(uint),1,fv);
		  fwrite(&(_repK2),sizeof(uint),1,fv);
		  fwrite(&(_maxRealLevel1),sizeof(uint),1,fv);
		  fwrite(&(_maxLevel1),sizeof(uint),1,fv);
		  fwrite(&(_maxLevel2),sizeof(uint),1,fv);
			printf("%d\n", _maxLevel2);

*/
		//  fwrite(&(trep->maxLevel),sizeof(uint),1,fv);
		//  
	fwrite(&zeroNode,sizeof(uint),1,fv); //stores the number of words of the vocabulary
	fwrite(&lenWords,sizeof(uint),1,fv);
	//Writes the vocabulary to disk.

	for (i=0;i<zeroNode;i++)
		fwrite(hash[positionInTH[i]].word,sizeof(char),lenWords,fv);

	fclose(fv);
	
	strcpy(filename,argv[1]);
	strcat(filename,".cil");
	FILE * fi = fopen(filename,"w");


	for(fila=0;fila<_part;fila++){
		for(columna=0;columna<_part;columna++){  
			  //fprintf(stderr,"fila: %d, columna: %d, edges: %d\n",fila,columna,_partnumberOfEdges[fila][columna]);
			fwrite(&(_partnumberOfNodes[fila][columna]),sizeof(uint),1,fi);
			fwrite(&(_partnumberOfEdges[fila][columna]),sizeof(ulong),1,fi);
			if(_partnumberOfEdges[fila][columna]==0)
				continue;
			fwrite(&(_partcutBt[fila][columna]),sizeof(uint),1,fi);
			fwrite(&(_partlastBt1_len[fila][columna]),sizeof(uint),1,fi);
			  //fwrite(&(rep->maxLevel),sizeof(uint),1,fi);
			fwrite(&(_partnleaves[fila][columna]),sizeof(uint),1,fi);
			  //fwrite (rep->leavesInf,sizeof(uint),_partnleaves[fila][columna]*K2*K2/W+1,fi);
			saveFT(_partcompressIL[fila][columna],fi);
		}
	}
	fclose(fi);   
	double t = 0;

	t += stop_clock(); 
  	t *= 1000; // to milliseconds
  	fprintf(stderr,"# compression_time = %f\n",t/1000);

	struct rusage resource_usage;
    getrusage(RUSAGE_SELF, &resource_usage);  		

  	fprintf(stderr,"# compression_space = %ld\n",resource_usage.ru_maxrss*1024);
  	fprintf(stderr,"# constructs_space_vmem = 0\n");
  	fprintf(stderr,"# construct_morton_duration = 0\n");
  	fprintf(stderr,"# construct_bv_complete_duration = 0\n");
  	fprintf(stderr,"# construct_sort_duration = 0\n");
  	fprintf(stderr,"# construct_duration = 0\n");
  	fprintf(stderr,"# buildvec_duration = 0\n");
  	fprintf(stderr,"# subtree_constructor_duration = 0\n");
  	fprintf(stderr,"# constructor_call_duration = 0\n");

	return 0;
}


