/*
    k-ary tree
        creation
        traversal  
        
    Notes: 
    1.  both creation and traversal follow BFS      
    2.  parent pointer is stored just to show that the tree is created
        correctly, its not used in any way to aid in creation and 
        traversal of the tree
*/

#include "kTree.h"

#define swap( x, y ) { int temp; temp=*(x); *(x)=*(y); *(y)=temp; }

typedef struct edgeinfo
{
	uint xval;
	uint yval;
	uint kval;
}tedge;


typedef struct QUEUE
{
	NODE *element;
	struct QUEUE *link;
	int basex;
	int basey;
}QUEUE;/*for implementing BFS*/

typedef struct QUEUE2
{
	int element;
	struct QUEUE2 *link;
	int basey;
	int basex;
}QUEUE2;/*for implementing BFS*/

typedef struct QUEUECONS
{
	NODE *element;
	struct QUEUECONS *link;

}QUEUECONS;/*for implementing BFS*/


typedef struct QUEUEOFFCONS
{
	ulong offsetL;
	ulong offsetR;
	struct QUEUEOFFCONS *link;
}QUEUEOFFCONS;

QUEUE * finalQUEUE;
QUEUE2 * finalQUEUE2;
QUEUECONS * finalQUEUECONS;
QUEUEOFFCONS * finalQUEUEOFFCONS;



uint *div_level_table;

unsigned int * positionInTH;
unsigned int addrInTH;
unsigned int zeroNode;

FILE * _ftr, * _flv,  * _fil;
ulong _ttlLeaves;

int hasAnyBitSet(unsigned char * b){
	int i=0,bvalue=0;
	if(b!=NULL)
		for(i=0;i<K2_3_char;i++)
			bvalue = bvalue||(b[i]!=0);
		return bvalue;
	}	


	QUEUE * AddItem (QUEUE * listpointer, NODE * elem, int cantx, int canty) {

  /*QUEUE * lp = listpointer;
  if (listpointer != NULL) {
		while (listpointer -> link != NULL)
	    listpointer = listpointer -> link;
		listpointer -> link = (QUEUE  *) malloc (sizeof (struct QUEUE));
		listpointer = listpointer -> link;
		listpointer -> link = NULL;
		listpointer -> element = elem;
		listpointer -> basex = cantx;
		listpointer -> basey = canty;
		return lp;
	}
  else {
		listpointer = (QUEUE  *) malloc (sizeof (struct QUEUE));
		listpointer -> link = NULL;
		listpointer -> element = elem;
		listpointer -> basex = cantx;
		listpointer -> basey = canty;
		return listpointer;
  }
  */
  if(listpointer != NULL){
  	QUEUE * lp = (QUEUE  *) malloc (sizeof (struct QUEUE));
  	finalQUEUE -> link = lp;
  	lp -> link = NULL;
  	lp -> element = elem;
  	lp -> basex = cantx;
  	lp -> basey = canty;
  	finalQUEUE = lp;
  	return listpointer;
  }
  else{
  	listpointer = (QUEUE  *) malloc (sizeof (struct QUEUE));
  	listpointer -> link = NULL;
  	listpointer -> element = elem;
  	listpointer -> basex = cantx;
  	listpointer -> basey = canty;
  	finalQUEUE = listpointer;
  	return listpointer;
  }

}

QUEUE * RemoveItem (QUEUE * listpointer) {

	QUEUE * tempp;
	tempp = listpointer -> link;
	free (listpointer);
	return tempp;
}

void ClearQueue (QUEUE * listpointer) {

	while (listpointer != NULL) {
		listpointer = (QUEUE *)RemoveItem (listpointer);
	}
}




void AddItem2 (TREP *trep, int elem, int cantx,int canty) {
	if(trep->iniq!=-1){
		trep->finq++;
		trep -> element[trep->finq] = elem;
		trep -> basex[trep->finq] = cantx;
		trep -> basey[trep->finq] = canty;

	}
	else{
		trep->iniq=0;
		trep->finq=0;
		trep -> element[trep->iniq] = elem;
		trep -> basex[trep->iniq] = cantx;
		trep -> basey[trep->iniq] = canty;
	}
}

void RemoveItem2 (TREP * trep) {

	trep->iniq++;
}

void AddItem3 (TREP *trep, int elem, int cantx,int canty, int p1, int p2, int q1, int q2) {
	if(trep->iniq!=-1){
		trep->finq++;
		trep -> element[trep->finq] = elem;
		trep -> basex[trep->finq] = cantx;
		trep -> basey[trep->finq] = canty;
		trep -> basep1[trep->finq]= p1;
		trep -> basep2[trep->finq]= p2;
		trep -> baseq1[trep->finq]= q1;
		trep -> baseq2[trep->finq]= q2;

	}
	else{
		trep->iniq=0;
		trep->finq=0;
		trep -> element[trep->iniq] = elem;
		trep -> basex[trep->iniq] = cantx;
		trep -> basey[trep->iniq] = canty;
		trep -> basep1[trep->finq]= p1;
		trep -> basep2[trep->finq]= p2;
		trep -> baseq1[trep->finq]= q1;
		trep -> baseq2[trep->finq]= q2;
	}
}

void RemoveItem3 (TREP * trep) {

	trep->iniq++;
}


QUEUECONS * AddItemCONS (QUEUECONS * listpointer, NODE * elem) {

	if(listpointer!=NULL){
		QUEUECONS * lp = (QUEUECONS  *) malloc (sizeof (struct QUEUECONS));
		finalQUEUECONS -> link = lp;
		lp -> link = NULL;
		lp -> element = elem;
		finalQUEUECONS = lp;
		return listpointer;
	}
	else{
		listpointer = (QUEUECONS  *) malloc (sizeof (struct QUEUECONS));
		listpointer -> link = NULL;
		listpointer -> element = elem;
		finalQUEUECONS = listpointer;
		return listpointer;
	}
}

QUEUECONS * RemoveItemCONS (QUEUECONS * listpointer) {
	QUEUECONS * tempp;
	tempp = listpointer -> link;
	free (listpointer);
	return tempp;
}


QUEUEOFFCONS * AddItemOFFCONS (QUEUEOFFCONS * listpointer, uint offsetL, uint offsetR) {

	if(listpointer!=NULL){
		QUEUEOFFCONS * lp = (QUEUEOFFCONS  *) malloc (sizeof (struct QUEUEOFFCONS));	
		finalQUEUEOFFCONS -> link = lp;
		lp -> link = NULL;
		lp -> offsetL = offsetL;
		lp -> offsetR = offsetR;
		finalQUEUEOFFCONS = lp;
			/*fprintf(stderr,"lp-oL oR: %d %d\n",lp -> offsetL,lp -> offsetR);
			fprintf(stderr,"lipo-oL oR: %d %d\n",listpointer -> offsetL,listpointer -> offsetR);
			fprintf(stderr,"lipo-oL oR: %d %d\n",listpointer -> link->offsetL,listpointer ->link-> offsetR);
			fprintf(stderr,"final-oL oR: %d %d\n",finalQUEUEOFFCONS -> offsetL,finalQUEUEOFFCONS -> offsetR);
			*/
			return listpointer;
		}
		else{
			listpointer = (QUEUEOFFCONS  *) malloc (sizeof (struct QUEUEOFFCONS));
			listpointer -> link = NULL;
			listpointer -> offsetL = offsetL;
			listpointer -> offsetR = offsetR;
			//fprintf(stderr,"lipo-oL oR: %d %d\n",listpointer -> offsetL,listpointer -> offsetR);

			finalQUEUEOFFCONS = listpointer;
			return listpointer;
		}
	}

	QUEUEOFFCONS * RemoveItemOFFCONS (QUEUEOFFCONS * listpointer) {
		QUEUEOFFCONS * tempp;
//    printf ("Element removed is %d\n", listpointer -> dataitem);
		tempp = listpointer -> link;
		if(listpointer)
			free (listpointer);
		return tempp;
	}



uint exp_pow(uint base, uint pow){
	uint i, result = 1;
	for(i=0;i<pow;i++)
		result*=base;

	return result;
}


/*------------------------------------------------------------------
 Function used by qsort to compare two elements a and b
------------------------------------------------------------------*/
 int comparaFrecListaDesc(const void *a, const void *b)
 {
 	unsigned int *left,*right;
 	left =  (unsigned int *) a;
 	right = (unsigned int *) b;
 	if ( hash[*left].weight < hash[*right].weight)
 		return 1;
 	if ( hash[*left].weight > hash[*right].weight)
 		return -1;

 	return 0;
 }

/*------------------------------------------------------------------
 Initilizes the structures used.
 ---------------------------------------------------------------- */
 void initialize (unsigned long sizeVoc){
	//unsigned long i;

 	_memMgr = createMemoryManager();
 	initialize_hash(sizeVoc);

 	positionInTH = (unsigned int*) malloc (sizeVoc * sizeof(unsigned int));
 	zeroNode = 0;
 }

 TREP * compactCreateKTree(char * filein, char * fileout, uint subm,int _K1, int _K2, int levelcut){
 	FILE *f;

 	uint nodes;
 	ulong edges;

 	ulong i,k,j,queuecont, conttmp;
	long z;
 	uint node,div_level, fila, columna;

 	uint tempk,tempx,tempy;

 	uint numberOfNodes;
 	ulong numberOfEdges;


 	f = fopen(filein,"r");
 	fread(&nodes,sizeof(uint),1,f);

 	uint part;
 	part = nodes/subm+1;


 	uint nodesOrig = nodes;
 	nodes= subm; 


 	uint max_level1;
 	uint max_real_level1;
 	max_real_level1 = ceil(log(nodes)/log(_K1))-1;
 	if(levelcut==0)
 		max_level1 = ceil(log(nodes)/log(_K1))-L;
 	else
 		max_level1 = levelcut;




 	uint nodes2 = 1;
 	for(i=0;i<max_real_level1+1;i++){
 		nodes2 *=_K1;
 	}


 	for(i=0;i<max_level1;i++)
 		nodes2 = ceil((double)nodes2/_K1);


 	uint	max_level2 = ceil(log(nodes2)/log(_K2))-1;  

 	uint maxreallevel1= max_real_level1;
 	uint maxlevel1 = max_level1;
 	uint maxlevel2 = max_level2;


 	fread(&edges,sizeof(ulong),1,f);

 	uint nodes_read=0;

 	uint algo;
 	ulong algo2;
 	int id1, id2;

 	uint kval, xedge, yedge;



 	ulong * counterP = (ulong *)malloc(sizeof(ulong)*part*part);

 	ulong * startP = (ulong *)malloc(sizeof(ulong)*part*part);
 	ulong * finP = (ulong *)malloc(sizeof(ulong)*part*part);

 	uint * nodesreadP = (uint *)malloc(sizeof(uint)*part*part);

 	for(j=0;j<part*part;j++){
 		counterP[j]=0;
 		startP[j]=0;
 		finP[j]=0;
 	}
 	counterP[0]=0;				


 	ulong edges_read=0;
 	uint cedg = 0;
 	ulong ultStart = 0;

 	uint * totaledges = (uint *)malloc(sizeof(uint)*(nodesOrig+edges));
	
 	for(i=0;i<nodesOrig+edges;i++) {
 		int k;
 		fread(&k,sizeof(uint),1,f);
 		totaledges[i]=k;
 		if(k<0) {
 			nodes_read++;
 			ultStart=i;
 		}
 		else {
 			k--;

 			xedge=nodes_read-1;
 			yedge=k;

 			kval = xedge/subm*part+yedge/subm;


 			if(counterP[kval]==0){
 				nodesreadP[kval]=nodes_read-1;
 				startP[kval]=ultStart;
 			}
 			finP[kval]=i;

 			counterP[kval]++;


 			cedg++;			     	
 			edges_read++;
 		}
 	}

 	fclose(f);


	
 	numberOfNodes = nodesOrig;
 	numberOfEdges = edges;
 	TREP * trep;
 	trep = malloc(sizeof(struct treeRep));
 	trep->part=part;
 	trep->tamSubm=subm;
 	trep->numberOfNodes=numberOfNodes;
 	trep->numberOfEdges=numberOfEdges;
 	trep->maxRealLevel1=maxreallevel1;
 	trep->maxLevel1=maxlevel1;
 	trep->maxLevel2=maxlevel2-L;
 	trep->repK1=_K1;
 	trep->repK2=_K2;

 	openPartialFile(trep,fileout);

/*	fwrite(&(trep->part),sizeof(uint),1,_fil);
 	fwrite(&(trep->tamSubm),sizeof(uint),1,_fil);

 	fwrite(&(trep->numberOfNodes),sizeof(uint),1,_fil);
 	fwrite(&(trep->numberOfEdges),sizeof(ulong),1,_fil);

 	fwrite(&(trep->repK1),sizeof(uint),1,_fil);
 	fwrite(&(trep->repK2),sizeof(uint),1,_fil);
 	fwrite(&(trep->maxRealLevel1),sizeof(uint),1,_fil);
 	fwrite(&(trep->maxLevel1),sizeof(uint),1,_fil);
 	fwrite(&(trep->maxLevel2),sizeof(uint),1,_fil);
 	*/




 	ulong nedges = numberOfEdges;
 	trep->submatrices=(MREP ***)malloc(sizeof(MREP **)*part);

 	for(i=0;i<part;i++){
 		trep->submatrices[i]=(MREP **)malloc(sizeof(MREP *)*part);
 	}



 	trep->info = (uint *)malloc(sizeof(uint)*MAX_INFO);			
 	trep->info2[0] = (uint *)malloc(sizeof(uint)*MAX_INFO);
 	trep->info2[1] = (uint *)malloc(sizeof(uint)*MAX_INFO);
 	trep->element = (uint *)malloc(sizeof(uint)*MAX_INFO);	
 	trep->basex = (uint *)malloc(sizeof(uint)*MAX_INFO);
 	trep->basey = (uint *)malloc(sizeof(uint)*MAX_INFO);

 	trep->iniq = -1;
 	trep->finq =-1;

 	trep->div_level_table1 = (uint *)malloc(sizeof(uint)*trep->maxLevel1);
 	for(i=0;i<trep->maxLevel1;i++){
 		trep->div_level_table1[i]=exp_pow(K1,trep->maxRealLevel1-i);
 	}

 	
 	trep->div_level_table2 = (uint *)malloc(sizeof(uint)*trep->maxLevel2);
 	for(i=0;i<trep->maxLevel2;i++){
 		trep->div_level_table2[i]=exp_pow(K2,trep->maxLevel2+L-i);
 	}

 	uint p;

 	_ttlLeaves=0;








 	for(fila=0;fila<part;fila++){
 		for(columna=0;columna<part;columna++){
 			ulong edges_sub=0;

 			uint p = fila*part+columna;

 			numberNodes=0;
 			numberLeaves=0;
 			numberTotalLeaves=0;

 			MREP * rep = trep->submatrices[fila][columna];
 			rep = malloc(sizeof(struct matrixRep));
 			rep->numberOfNodes = subm;

 			rep->numberOfEdges = counterP[fila*part+columna];

 			nedges = rep->numberOfEdges;


//			fprintf(stderr,"creating k2-tree... fila: %d columna: %d part: %d nedges: %lld\n",fila,columna,part,nedges);
			fprintf(stderr,"creating k2-tree... %2.2f nodes processed\n",(float)(fila*part+columna)*100/part/part);
//			fprintf(stderr,"row= %d, column=%d, edges from %d and from %d, (%2.2f\%)\n",fila,columna,fila*subm, columna*subm,(float)(fila*part+columna)*100/part/part);


 			nodes_read=nodesreadP[p];
 			edges_read=0;

 			tedge * tedges =(tedge *) malloc(sizeof(tedge)*nedges);
			if(nedges>0){

 			for(i=startP[p];i<=finP[p];i++) {
 				int k=totaledges[i];
 				if(k<0) {
 					nodes_read++;
 				}
 				else {
 					k--;


 					id1=nodes_read-1;
 					id2=k;

 					if((id1/subm ==fila)&&(id2/subm==columna)){

 						tedges[edges_sub].xval=id1%subm;
 						tedges[edges_sub].yval=id2%subm;
 						edges_sub++;
 					}
 					edges_read++;
 				}
 			}
//fprintf(stderr,"fila %d columna: %d p %d edges_read %lld counter %lld\n",fila,columna,p,edges_sub,counterP[p]);
 			if(edges_sub!=counterP[p]){
 				//fprintf(stderr,"fila %d columna: %d p %d edges_read %d counter %d\n",fila,columna,p,edges_sub,counterP[p]);
 				exit(0);
 			}


 			ulong sizebt = 0;
 			ulong sizebn = 0;
 			ulong sizeli = 0;
 			uint l = 0;

		ulong offsetL=0;//boundariesP[fila*part+columna]; 
		ulong offsetR=counterP[p];//boundariesP[fila*part+columna+1];

		ulong postotal=0;

		QUEUEOFFCONS * q = NULL;
		finalQUEUEOFFCONS = q;
		
		uint bits_BT_len = 0;
		uint bits_BN_len = 0;
		uint bits_LI_len = 0;
		
		ulong * counterK1 = (ulong *)malloc(sizeof(ulong)*K1*K1);
		ulong * boundariesK1 = (ulong *)malloc(sizeof(ulong)*(K1*K1+1));
		ulong * pointerK1 = (ulong *)malloc(sizeof(ulong)*(K1*K1));


		ulong * counterK2 = (ulong *)malloc(sizeof(ulong)*K2*K2);
		ulong * boundariesK2 = (ulong *)malloc(sizeof(ulong)*(K2*K2+1));
		ulong * pointerK2 = (ulong *)malloc(sizeof(ulong)*(K2*K2));


		ulong * counterK2_2 = (ulong *)malloc(sizeof(ulong)*K2_2*K2_2);



		uint tempk,tempx,tempy;


		ulong lenBT;
		uint * bits;
		bitRankW32Int *BT, *BN;

		uint * bits_BN;
		uint tamBT = nedges;
		if(nedges<nodesreadP[p])
			tamBT = nodesreadP[p];
		uint * bits_BT=(uint *) malloc(sizeof(uint)*tamBT);
		for(i=0;i<tamBT;i++)
			bits_BT[i]=0;
		bits = bits_BT;

		//fprintf(stderr,"nedges: %d\n",nedges);	
		q = AddItemOFFCONS(q,offsetL,offsetR);
		queuecont = 1;
		
		
		uint K= trep->repK1;

		ulong * counterK = counterK1;
		ulong * pointerK = pointerK1;
		ulong * boundariesK = boundariesK1;

		
		for(i=0;i<trep->maxLevel1+trep->maxLevel2;i++){
			conttmp = 0;
			if(i<trep->maxLevel1)
				div_level = trep->div_level_table1[i];
			
			if(i==trep->maxLevel1-1) {					
				rep->lastBt1_len = postotal;
			}
			
			if(i==trep->maxLevel1) {					
				rep->cutBt = postotal;
			//	fprintf(stderr,"pos para cutBt: %d lastBt1_len %d\n",rep->cutBt,rep->lastBt1_len);
			}
			
			
			if(i>=trep->maxLevel1){
				K=trep->repK2;
				counterK = counterK2;
				pointerK = pointerK2;
				boundariesK = boundariesK2;

				div_level = trep->div_level_table2[i-trep->maxLevel1];

			}
			if(i==trep->maxLevel1+trep->maxLevel2-1){
				bits_BT_len=postotal;
				postotal=0;
				bits_BN_len= queuecont*K*K;
				bits_BN = (uint*)malloc(sizeof(uint)*((bits_BN_len/W+1)));
				for(j=0; j<(bits_BN_len/W+1);j++)
					bits_BN[j]=0;
				bits = bits_BN;

			}
			
			for(k=0;k<queuecont;k++){
				
				if(i<trep->maxLevel1+trep->maxLevel2-1)
					sizebt +=K*K;
				else 
					sizebn +=K*K;
				
				
				offsetL = q->offsetL;
				offsetR = q->offsetR;
				
				for(j=0;j<K*K;j++){
					counterK[j]=0;
					pointerK[j]=0;
				}

				for(j=offsetL;j<offsetR;j++){
					if(offsetL>offsetR){
						fprintf(stderr,"j: %d divlevel: %d\n",j, div_level);
						fprintf(stderr,"%d %d\n", tedges[j].xval , tedges[j].yval );
						}

					tedges[j].kval= (tedges[j].xval / div_level)*K+tedges[j].yval /div_level;
					tedges[j].xval=tedges[j].xval%div_level;
					tedges[j].yval=tedges[j].yval%div_level;
					if(tedges[j].kval>=K*K){
						fprintf(stderr,"ERROR ");
						fprintf(stderr,"K %d  xval %d yval %d kval %d\n",K,tedges[j].xval,tedges[j].yval,tedges[j].kval);
						}
					counterK[tedges[j].kval]++;
				}
				
				boundariesK[0]=offsetL;			
				for(j=0;j<K*K;j++){
					boundariesK[j+1]=boundariesK[j]+counterK[j];
					pointerK[j]=boundariesK[j];

					if(boundariesK[j+1]!=boundariesK[j]){
						conttmp++;
						q=AddItemOFFCONS(q,boundariesK[j],boundariesK[j+1]);
						bitset(bits,postotal);
					}
					postotal++;

				}

				for(z=0;z<K*K;z++){
					while(pointerK[z]<boundariesK[z+1]){
						if(tedges[pointerK[z]].kval!=z){
							tempk = tedges[pointerK[z]].kval;
							tempx = tedges[pointerK[z]].xval;
							tempy = tedges[pointerK[z]].yval;

							while(tedges[pointerK[tempk]].kval==tempk){
								pointerK[tempk]++;
							}

							tedges[pointerK[z]].kval = tedges[pointerK[tempk]].kval;
							tedges[pointerK[z]].xval = tedges[pointerK[tempk]].xval;
							tedges[pointerK[z]].yval = tedges[pointerK[tempk]].yval;

							tedges[pointerK[tempk]].kval= tempk;
							tedges[pointerK[tempk]].xval= tempx;
							tedges[pointerK[tempk]].yval= tempy;
							pointerK[tempk]++;


						}
						else
							pointerK[z]++;
						
					}
				}

				
				for(z=offsetL;z<(long)offsetR-1;z++)
					if(tedges[z].kval>tedges[z+1].kval){
						fprintf(stderr,"Error? z %d %d %d offsets: %lld %lld, i %d, p: %d\n",z,tedges[z].kval,tedges[z+1].kval,offsetL, offsetR,i,p);
						exit(0);
					}
					q = (QUEUEOFFCONS *)RemoveItemOFFCONS(q);
				}
				queuecont = conttmp;
			}

			K = K2_2;
			counterK = counterK2_2;
			bits_LI_len = queuecont*K*K;
			uint * bits_LI = (uint*)malloc(sizeof(uint)*((bits_LI_len+W-1)/W));
			for(i=0; i<(W-1+bits_LI_len)/W;i++)
				bits_LI[i]=0;

			uint counttotal=0;
			postotal=0;
			while(q!=NULL){
				sizeli +=K*K;
				
				offsetL = q->offsetL;
				offsetR = q->offsetR;

				for(j=0;j<K*K;j++){
					counterK[j]=0;
				}
				
				div_level = K;	
				for(j=offsetL;j<offsetR;j++){
					tedges[j].kval= tedges[j].xval *K+tedges[j].yval;
					counterK[tedges[j].kval]++;
				}

				
				for(j=0;j<K*K;j++){

					if(counterK[j]>0){
						conttmp++;
						counttotal++;

						bitset(bits_LI,postotal);
					}
					postotal++;
					
					if(counterK[j]>1){
						fprintf(stderr,"ofL oR x y: %d %d %d %d k: %d\n",offsetL, offsetR, tedges[j].xval,tedges[j].yval,tedges[j].kval );
						fprintf(stderr,"error : j %d counter: %d \n",j,counterK[j]);
					}

				}

				q = (QUEUEOFFCONS *)RemoveItemOFFCONS(q);
			}



			free(counterK);
			free(boundariesK);
			free(pointerK);

			BT = createBitRankW32Int(bits_BT, bits_BT_len , 1, 20); 


			rep->bt = BT;
			rep->bt_len = bits_BT_len;


			BN = createBitRankW32Int(bits_BN, bits_BN_len , 1, 20); 


			rep->bn = BN;
			rep->bn_len = bits_BN_len;

			rep->leavesInf = bits_LI;
			rep->nleaves = bits_LI_len/(K*K);
			
		//PARTIAL SAVE

			if(rep->numberOfEdges!=0)
				_ttlLeaves +=rep->nleaves;	


			if(rep->numberOfEdges!=0)
				save(rep->bt,_ftr);


			if(rep->numberOfEdges!=0)
				save(rep->bn,_flv);
			
//			fprintf(stderr,"nodes: %d edges: %d\n",rep->numberOfNodes,rep->numberOfEdges);

			destroyBitRankW32Int(rep->bt);
			destroyBitRankW32Int(rep->bn);
			}
			fwrite(&(rep->numberOfNodes),sizeof(uint),1,_fil);
			fwrite(&(rep->numberOfEdges),sizeof(ulong),1,_fil);
			if(rep->numberOfEdges>0){

				fwrite(&(rep->cutBt),sizeof(uint),1,_fil);
				fwrite(&(rep->lastBt1_len),sizeof(uint),1,_fil);
				fwrite(&(rep->nleaves),sizeof(uint),1,_fil);
				if(rep->numberOfEdges==0)
					return;
				fwrite(rep->leavesInf,sizeof(uint),(rep->nleaves*K2_2*K2_2+W-1)/W,_fil);
				free(rep->leavesInf);
			}


			trep->submatrices[fila][columna]=rep;

			free(tedges);
			//fprintf(stderr,"fin de fila %d columna %d\n",fila,columna);

		}
	}

	free(nodesreadP);
	free(startP);
	free(finP);
	free(counterP);
	free(totaledges);
	return trep;
}


NODE * createKTree(int _K1, int _K2, int maxreallevel1, int maxlevel1, int maxlevel2){
	NODE * n = (NODE *) malloc(sizeof(struct node));
	n->child=NULL;
	n->data=0;
//	K1 = _K1;
//	K2 = _K2;
	max_real_level1 = maxreallevel1;
	max_Level1 = maxlevel1;
	max_Level2 = maxlevel2-L;
	numberNodes =0;
	numberLeaves = 0;
	numberTotalLeaves=0;
	return n;
}


MREP * loadRepresentation(char * basename){
	MREP * rep;
	int i;
	rep = (MREP *) malloc(sizeof(struct matrixRep));
	rep->bt = (bitRankW32Int *) malloc(sizeof(struct sbitRankW32Int));
	rep->bn = (bitRankW32Int *) malloc(sizeof(struct sbitRankW32Int));	

	char *filename = (char *) malloc(sizeof(char)*(strlen(basename)+4));
	strcpy(filename,basename);
	strcat(filename,".tr");
	FILE * ft = fopen(filename,"r");
	load(rep->bt,ft);
	fclose(ft);
	rep->bt_len = rep->bt->n;

	strcpy(filename,basename);
	strcat(filename,".lv");
	FILE * fn = fopen(filename,"r");
	load(rep->bn,fn);
	fclose(fn);  
	rep->bn_len = rep->bn->n;



	strcpy(filename,basename);
	strcat(filename,".il");
	FILE * fi = fopen(filename,"r");
	fread(&(rep->numberOfNodes),sizeof(uint),1,fi);

	fread(&(rep->numberOfEdges),sizeof(ulong),1,fi);


	fread(&(rep->cutBt),sizeof(uint),1,fi);

	fread(&(rep->lastBt1_len),sizeof(uint),1,fi);

	fread(&(rep->nleaves),sizeof(uint),1,fi);

	rep->leavesInf = (uint *)malloc(sizeof(uint)*(rep->nleaves*K2*K2/W+1));
	fread(rep->leavesInf,sizeof(uint),rep->nleaves*K2*K2/W+1,fi);
	fclose(fi);   


	free(filename);
	return rep;
}

void destroyRepresentation(MREP * rep){
	destroyBitRankW32Int(rep->bt);
	destroyBitRankW32Int(rep->bn);
	destroyFT(rep->compressIL);
	free(rep);
}



uint * compactAdjacencyList(TREP * trep,MREP * rep, int x){
	int K1K1 = K1*K1;
	int K2K2 = K2*K2;
	
	if(rep->numberOfEdges==0)
		return trep->info;

	uint factorAdjust = K1K1/K2K2;
	trep->iniq=-1;
	trep->finq=-1;
	uint nleaf,posInf,realvalue,  nleafrelat;

	uint totalAdyNodes =0;
	int i, k, j, queuecont, conttmp,node,div_level,xrelat;

	AddItem2(trep,0,0,x);

	queuecont = 1;
	for(i=0;i<trep->maxLevel1-1;i++){
		conttmp = 0;


		div_level = trep->div_level_table1[i];

		for(k=0;k<queuecont;k++){
			for(j=0;j<K1;j++){
				xrelat = trep->basey[trep->iniq];
				node = xrelat/div_level*K1 + j;
				node += trep->element[trep->iniq];


				if(isBitSet(rep->bt,node)){

					conttmp++;

					AddItem2(trep,rank(rep->bt,node)*K1K1,trep->basex[trep->iniq]+j*div_level,trep->basey[trep->iniq]%div_level);

				}
			}
			
			RemoveItem2(trep);
		}
		queuecont = conttmp;
	}


	for(i=trep->maxLevel1-1;i<trep->maxLevel1;i++){
		conttmp = 0;


		div_level = trep->div_level_table1[i];

		for(k=0;k<queuecont;k++){
			for(j=0;j<K1;j++){
				xrelat = trep->basey[trep->iniq];
				node = xrelat/div_level*K1 + j;
				node += trep->element[trep->iniq];


				if(isBitSet(rep->bt,node)){

					conttmp++;

					AddItem2(trep,rep->cutBt + (rank(rep->bt,node)*K1K1-rep->cutBt)/factorAdjust,trep->basex[trep->iniq]+j*div_level,trep->basey[trep->iniq]%div_level);

				}
			}
			
			RemoveItem2(trep);
		}
		queuecont = conttmp;
	}

	uint cutPreRank = rank(rep->bt,rep->lastBt1_len-1);

	for(i=0;i<trep->maxLevel2-1;i++){
		conttmp = 0;


		div_level = trep->div_level_table2[i];

		for(k=0;k<queuecont;k++){
			for(j=0;j<K2;j++){
				xrelat = trep->basey[trep->iniq];
				node = xrelat/div_level*K2 + j;
				node += trep->element[trep->iniq];


				if(isBitSet(rep->bt,node)){

					conttmp++;

					AddItem2(trep,rep->cutBt + (rank(rep->bt,node-1)-cutPreRank)*K2K2,trep->basex[trep->iniq]+j*div_level,trep->basey[trep->iniq]%div_level);

				}
			}
			
			RemoveItem2(trep);
		}
		queuecont = conttmp;
	}

	
	while(trep->iniq<=trep->finq){
		nleaf = trep->element[trep->iniq]-rep->bt_len;

		for(j=0;j<K2;j++){
			nleafrelat = nleaf + (trep->basey[trep->iniq]/K2_2)*K2+j;

			if(isBitSet(rep->bn,nleafrelat)){

				posInf = rank(rep->bn,nleafrelat);

				realvalue = accessFT(rep->compressIL,posInf);				
				for(i=0;i<K2_2;i++){

					if(bitgetchar(&(trep->words[realvalue*trep->lenWords]),(i+(x%K2_2)*K2_2))){

						trep->info[0]++;
						trep->info[trep->info[0]]=trep->basex[trep->iniq]+i+K2_2*j+trep->columna*trep->tamSubm;
					}
				}
			}
		}
		RemoveItem2(trep);
	}

	return trep->info;
}


uint * compactInverseList(TREP * trep,MREP * rep, int y){
	if(rep->numberOfEdges==0)
		return trep->info;

	int K1K1 = K1*K1;
	int K2K2 = K2*K2;
	uint factorAdjust = K1K1/K2K2;
	trep->iniq=-1;
	trep->finq=-1;
	uint nleaf,posInf,realvalue,  nleafrelat;

	uint totalAdyNodes =0;
	int i, k, j, queuecont, conttmp,node,div_level,yrelat;

	AddItem2(trep,0,y,0);

	queuecont = 1;
	for(i=0;i<trep->maxLevel1-1;i++){
		conttmp = 0;


		div_level = trep->div_level_table1[i];

		for(k=0;k<queuecont;k++){
			for(j=0;j<K1;j++){
				yrelat = trep->basex[trep->iniq];
				node = K1*j + yrelat/div_level;
				node += trep->element[trep->iniq];


				if(isBitSet(rep->bt,node)){

					conttmp++;

					AddItem2(trep,rank(rep->bt,node)*K1K1,trep->basex[trep->iniq]%div_level,trep->basey[trep->iniq]+j*div_level);

				}
			}
			
			RemoveItem2(trep);
		}
		queuecont = conttmp;
	}


	for(i=trep->maxLevel1-1;i<trep->maxLevel1;i++){
		conttmp = 0;


		div_level = trep->div_level_table1[i];

		for(k=0;k<queuecont;k++){
			for(j=0;j<K1;j++){
				yrelat = trep->basex[trep->iniq];
				node = K1*j+ yrelat/div_level;
				node += trep->element[trep->iniq];


				if(isBitSet(rep->bt,node)){

					conttmp++;

					AddItem2(trep,rep->cutBt + (rank(rep->bt,node)*K1K1-rep->cutBt)/factorAdjust,trep->basex[trep->iniq]%div_level,trep->basey[trep->iniq]+j*div_level);

				}
			}
			
			RemoveItem2(trep);
		}
		queuecont = conttmp;
	}

	uint cutPreRank = rank(rep->bt,rep->lastBt1_len-1);

	for(i=0;i<trep->maxLevel2-1;i++){
		conttmp = 0;

		div_level = trep->div_level_table2[i];

		for(k=0;k<queuecont;k++){
			for(j=0;j<K2;j++){
				yrelat = trep->basex[trep->iniq];
				node = K2*j + yrelat/div_level;
				node += trep->element[trep->iniq];

				if(isBitSet(rep->bt,node)){
					conttmp++;
					AddItem2(trep,rep->cutBt + (rank(rep->bt,node-1)-cutPreRank)*K2K2,trep->basex[trep->iniq]%div_level,trep->basey[trep->iniq]+j*div_level);
				}
			}
			
			RemoveItem2(trep);
		}
		queuecont = conttmp;
	}

	
	while(trep->iniq<=trep->finq){
		nleaf = trep->element[trep->iniq]-rep->bt_len;
		for(j=0;j<K2;j++){
			nleafrelat = nleaf + (trep->basex[trep->iniq]/K2_2)+K2*j;
			if(isBitSet(rep->bn,nleafrelat)){
				posInf = rank(rep->bn,nleafrelat);
				realvalue = accessFT(rep->compressIL,posInf);

				for(i=0;i<K2_2;i++){
					if(bitgetchar(&(trep->words[realvalue*trep->lenWords]),(i*K2_2+(y%K2_2)))){
						trep->info[0]++;
						trep->info[trep->info[0]]=trep->basey[trep->iniq]+i+K2_2*j+trep->fila*trep->tamSubm;
					}
				}
			}
		}
		RemoveItem2(trep);
	}

	return trep->info;
}




MREP * createRepresentation(NODE * root, uint numberOfNodes,ulong numberOfEdges){
	MREP * rep;
	rep = malloc(sizeof(struct matrixRep));
	rep->numberOfNodes = numberOfNodes;
	rep->numberOfEdges = numberOfEdges;
	uint bits_BT_len = numberNodes;
	uint bits_BN_len = numberTotalLeaves;
	uint bits_LI_len = numberLeaves*K2_2*K2_2;
	bitRankW32Int *BT, *BN;
	uint * bits_BT = (uint*)malloc(sizeof(uint)*((bits_BT_len+W-1)/W));
	uint * bits_BN = (uint*)malloc(sizeof(uint)*((bits_BN_len+W-1)/W));
	uint * bits_LI = (uint*)malloc(sizeof(uint)*((bits_LI_len+W-1)/W));
	

	
	int i, k, j, queuecont, conttmp,node,div_level, pos=0;
	for(i=0; i<(W-1+bits_BT_len)/W;i++)
		bits_BT[i]=0;
	for(i=0; i<(W-1+bits_BN_len)/W;i++)
		bits_BN[i]=0;
	for(i=0; i<(W-1+bits_LI_len)/W;i++)
		bits_LI[i]=0;

	char isroot=1;
	QUEUECONS * q=NULL;
	finalQUEUECONS = q;
	q = AddItemCONS(q,root);
	queuecont = 1;
	for(i=0;i<max_Level1;i++){
		conttmp = 0;
		div_level = exp_pow(K1,max_real_level1-i);
		for(k=0;k<queuecont;k++){
			if(q->element->child!=NULL){
				for(j=0;j<K1*K1;j++){
					node = j;
					conttmp++;
					q=AddItemCONS(q,q->element->child[node]);
					
				}
				if(!isroot)
					bitset(bits_BT,pos);
				free(q->element->child);

			}
			if(!isroot)
				pos++;
			isroot=0;
			free(q->element);
			q = (QUEUECONS *)RemoveItemCONS(q);
		}
		queuecont = conttmp;
	}
	rep->lastBt1_len = pos;
	rep->cutBt = pos+queuecont;
	for(i=0;i<max_Level2;i++){
		conttmp = 0;
		div_level = exp_pow(K2,max_Level2+L-i);
		for(k=0;k<queuecont;k++){
			if(q->element->child!=NULL){
				for(j=0;j<K2*K2;j++){
					node = j;
					conttmp++;
					q=AddItemCONS(q,q->element->child[node]);
					
				}
				if(!isroot)
					bitset(bits_BT,pos);
				free(q->element->child);

			}
			pos++;
			free(q->element);			
			q = (QUEUECONS *)RemoveItemCONS(q);
		}
		queuecont = conttmp;
	}


	BT = createBitRankW32Int(bits_BT, bits_BT_len , 1, 20); 


	rep->bt = BT;
	rep->bt_len = pos;

	pos=0;
	uint pos_inf=0;
	while(q!=NULL){
		if(hasAnyBitSet((q->element)->data)){
			bitset(bits_BN,pos);
			
			for(i=0;i<K2_2*K2_2;i++){
				if(bitgetchar((q->element)->data,i)){
					bitset(bits_LI,pos_inf);
				}
				pos_inf++;
			}
		}
		pos++;
		free(q->element->data);
		free(q->element);
		q = (QUEUECONS *)RemoveItemCONS(q);
	}

	BN = createBitRankW32Int(bits_BN, bits_BN_len , 1, 20); 

	rep->bn = BN;
	rep->bn_len = pos;
	
	rep->leavesInf = bits_LI;
	rep->nleaves = numberLeaves;
	
	return rep;
}

void insertNode(NODE * root, int x, int y){
	uint i,node, div_level;
	int l=0;
	NODE * n = root;
	while(l<max_Level1){
		div_level = exp_pow(K1,max_real_level1-l);
		node = (x / div_level)*K1+y /div_level;
		if(n->child==NULL){
			numberNodes+=K1*K1;
			n->child = (NODE **)malloc(sizeof(NODE *)*K1*K1);
			for(i=0;i<K1*K1;i++){
				n->child[i]=(NODE *) malloc(sizeof(struct node));
				n->child[i]->child=NULL;
				n->child[i]->data=0;
			}
		}
		n = n->child[node];

		
		x = x % div_level;
		y = y % div_level;
		l++;
	}
	
	l=0;
	int j;
	while(l<=max_Level2){
		div_level = exp_pow(K2,max_Level2+L-l);
		node = (x / div_level)*K2+y /div_level;
		if(l==max_Level2){
			if(!hasAnyBitSet(n->data)){
				numberLeaves++;
			}
			node = x *K2_2+y ;
			bitsetchar(n->data,node);
		}
		else{
			if(n->child==NULL){
				if(l<max_Level2-1)
					numberNodes+=K2*K2;
				else
					numberTotalLeaves+=K2*K2;
				n->child = (NODE **)malloc(sizeof(NODE *)*K2*K2);
				for(i=0;i<K2*K2;i++){
					n->child[i]=(NODE *) malloc(sizeof(struct node));
					n->child[i]->child=NULL;
					
					if(l==max_Level2-1){
						n->child[i]->data=malloc(sizeof(char)*K2_3_char);
						for(j=0;j<K2_3_char;j++)
							n->child[i]->data[j]=0;
					}
					
				}
			}
			n = n->child[node];
			
		} 
		
		x = x % div_level;
		y = y % div_level;
		l++;
	}

	
}




TREP * createTreeRep(uint nodesOrig,ulong edges,uint part,uint subm, uint max_real_level1, uint max_level1, uint max_level2, uint _K1, uint _K2){
	TREP * trep;
	trep = malloc(sizeof(struct treeRep));
	trep->part=part;
	trep->tamSubm=subm;
	trep->numberOfNodes=nodesOrig;
	trep->numberOfEdges=edges;
    trep->maxRealLevel1=max_real_level1;//-L;
    trep->maxLevel1=max_level1;
    trep->maxLevel2=max_level2-L;
    trep->repK1=_K1;
    trep->repK2=_K2;
    
    trep->submatrices=(MREP ***)malloc(sizeof(MREP **)*part);
    int i,j;
    for(i=0;i<part;i++){
    	trep->submatrices[i]=(MREP **)malloc(sizeof(MREP *)*part);
    }



    trep->info = (uint *)malloc(sizeof(uint)*MAX_INFO);

    trep->info2[0] = (uint *)malloc(sizeof(uint)*MAX_INFO);

    trep->info2[1] = (uint *)malloc(sizeof(uint)*MAX_INFO);

    trep->element = (uint *)malloc(sizeof(uint)*MAX_INFO);	

    trep->basex = (uint *)malloc(sizeof(uint)*MAX_INFO);

    trep->basey = (uint *)malloc(sizeof(uint)*MAX_INFO);

    trep->iniq = -1;
    trep->finq =-1;
    trep->div_level_table1 = (uint *)malloc(sizeof(uint)*trep->maxLevel1);
    for(i=0;i<trep->maxLevel1;i++){
    	trep->div_level_table1[i]=exp_pow(K1,trep->maxRealLevel1-i);
    }

    trep->div_level_table2 = (uint *)malloc(sizeof(uint)*trep->maxLevel2);
    for(i=0;i<trep->maxLevel2;i++){
    	trep->div_level_table2[i]=exp_pow(K2,trep->maxLevel2+L-i);

    }

    return trep;

}

void insertIntoTreeRep(TREP * trep, MREP * rep, uint i, uint j){
	trep->submatrices[i][j]=rep;
	
}


void saveTreeRep(TREP * trep, char * basename){
	char *filename = (char *) malloc(sizeof(char)*(strlen(basename)+5));
	strcpy(filename,basename);
	strcat(filename,".tr");
	FILE * ft = fopen(filename,"w");
	uint part= trep->part;
	int fila,columna;
	MREP* rep;
	for(fila=0;fila<part;fila++){
		for(columna=0;columna<part;columna++){
			if(trep->submatrices[fila][columna]->numberOfEdges!=0)
				save(trep->submatrices[fila][columna]->bt,ft);
		}
	}
	fclose(ft);


	strcpy(filename,basename);
	strcat(filename,".lv");
	FILE * fn = fopen(filename,"w");
	for(fila=0;fila<part;fila++){
		for(columna=0;columna<part;columna++){
			if(trep->submatrices[fila][columna]->numberOfEdges!=0)
				save(trep->submatrices[fila][columna]->bn,fn);
		}
	}
	fclose(fn);  


	strcpy(filename,basename);
	strcat(filename,".voc");
	FILE * fv = fopen(filename,"w");

	fwrite(&(trep->part),sizeof(uint),1,fv);
	fwrite(&(trep->tamSubm),sizeof(uint),1,fv);

	fwrite(&(trep->numberOfNodes),sizeof(uint),1,fv);
	fwrite(&(trep->numberOfEdges),sizeof(ulong),1,fv);

	fwrite(&(trep->repK1),sizeof(uint),1,fv);
	fwrite(&(trep->repK2),sizeof(uint),1,fv);
	fwrite(&(trep->maxRealLevel1),sizeof(uint),1,fv);
	fwrite(&(trep->maxLevel1),sizeof(uint),1,fv);
	fwrite(&(trep->maxLevel2),sizeof(uint),1,fv);


  fwrite(&trep->zeroNode,sizeof(uint),1,fv); //stores the number of words of the vocabulary
  fwrite(&trep->lenWords,sizeof(uint),1,fv);
	//Writes the vocabulary to disk.
  int i;
  for (i=0;i<zeroNode;i++)
  	fwrite(hash[positionInTH[i]].word,sizeof(char),trep->lenWords,fv);

  fclose(fv);

  strcpy(filename,basename);
  strcat(filename,".cil");
  FILE * fi = fopen(filename,"w");


  for(fila=0;fila<part;fila++){
  	for(columna=0;columna<part;columna++){  
  		rep=trep->submatrices[fila][columna];
  		fwrite(&(rep->numberOfNodes),sizeof(uint),1,fi);
  		fwrite(&(rep->numberOfEdges),sizeof(ulong),1,fi);
  		if(trep->submatrices[fila][columna]->numberOfEdges==0)
  			continue;
  		fwrite(&(rep->cutBt),sizeof(uint),1,fi);
  		fwrite(&(rep->lastBt1_len),sizeof(uint),1,fi);
  		fwrite(&(rep->nleaves),sizeof(uint),1,fi);
  		saveFT(rep->compressIL,fi);
  	}
  }
  fclose(fi);   
  
  

  
  
  free(filename);


}



TREP * loadTreeRepresentation(char * basename){
	TREP * trep;
	MREP * rep;
	int i,j,k;
	
	trep = malloc(sizeof(struct treeRep));
	char *filename = (char *) malloc(sizeof(char)*(strlen(basename)+5));


	strcpy(filename,basename);
	strcat(filename,".voc");
	FILE * fv = fopen(filename,"r");

	fread(&(trep->part),sizeof(uint),1,fv);
	fread(&(trep->tamSubm),sizeof(uint),1,fv);

	fread(&(trep->numberOfNodes),sizeof(uint),1,fv);
	fread(&(trep->numberOfEdges),sizeof(ulong),1,fv);
	
	
	fread(&(trep->repK1),sizeof(uint),1,fv);
	fread(&(trep->repK2),sizeof(uint),1,fv);
	fread(&(trep->maxRealLevel1),sizeof(uint),1,fv);
	fread(&(trep->maxLevel1),sizeof(uint),1,fv);
	fread(&(trep->maxLevel2),sizeof(uint),1,fv);

	
//	fread(&(trep->maxLevel),sizeof(uint),1,fv);
	
	fread(&trep->zeroNode,sizeof(uint),1,fv); //reads the number of words of the vocabulary
	fread(&trep->lenWords,sizeof(uint),1,fv);
	
	trep->words = (unsigned char *) malloc(sizeof(char)*trep->lenWords*trep->zeroNode);
	fread(trep->words,sizeof(unsigned char),trep->lenWords*trep->zeroNode,fv);

	
	fclose(fv);
	
	
	trep->div_level_table1 = (uint *)malloc(sizeof(uint)*trep->maxLevel1);
	for(i=0;i<trep->maxLevel1;i++)
		trep->div_level_table1[i]=exp_pow(K1,trep->maxRealLevel1-i);
	
	trep->div_level_table2 = (uint *)malloc(sizeof(uint)*trep->maxLevel2);
	for(i=0;i<trep->maxLevel2;i++)
		trep->div_level_table2[i]=exp_pow(K2,trep->maxLevel2+L-i);

	



	trep->info = (uint *)malloc(sizeof(uint)*MAX_INFO);
	trep->info2[0] = (uint *)malloc(sizeof(uint)*MAX_INFO);
	trep->info2[1] = (uint *)malloc(sizeof(uint)*MAX_INFO);
	trep->element = (uint *)malloc(sizeof(uint)*MAX_INFO);	
	trep->basex = (uint *)malloc(sizeof(uint)*MAX_INFO);
	trep->basey = (uint *)malloc(sizeof(uint)*MAX_INFO);
	trep->basep1 = (uint *)malloc(sizeof(uint)*MAX_INFO);
	trep->basep2 = (uint *)malloc(sizeof(uint)*MAX_INFO);
	trep->baseq1 = (uint *)malloc(sizeof(uint)*MAX_INFO);
	trep->baseq2 = (uint *)malloc(sizeof(uint)*MAX_INFO);

	trep->iniq = -1;
	trep->finq =-1;

	trep->submatrices=(MREP ***)malloc(sizeof(MREP **)*trep->part);

	int fila,columna;
	for(i=0;i<trep->part;i++){
		trep->submatrices[i]=(MREP **)malloc(sizeof(MREP *)*trep->part);
	}

	for(i=0;i<trep->part;i++){
		for(j=0;j<trep->part;j++){
			rep = (MREP *) malloc(sizeof(struct matrixRep));
			rep->bt = (bitRankW32Int *) malloc(sizeof(struct sbitRankW32Int));
			rep->bn = (bitRankW32Int *) malloc(sizeof(struct sbitRankW32Int));				
			trep->submatrices[i][j]=rep;
		}
	}

	strcpy(filename,basename);
	strcat(filename,".cil");
	FILE * fi = fopen(filename,"r");
	long sumatotalnleaves=0;
	for(fila=0;fila<trep->part;fila++){
		for(columna=0;columna<trep->part;columna++){  

			rep=trep->submatrices[fila][columna];
			fread(&(rep->numberOfNodes),sizeof(uint),1,fi);
			fread(&(rep->numberOfEdges),sizeof(ulong),1,fi);
			if(rep->numberOfEdges==0)
				continue;

			fread(&(rep->cutBt),sizeof(uint),1,fi);

			fread(&(rep->lastBt1_len),sizeof(uint),1,fi);


			fread(&(rep->nleaves),sizeof(uint),1,fi);

			sumatotalnleaves+=rep->nleaves;
			rep->compressIL = loadFT(fi);



		}
	}
	fclose(fi);   


	strcpy(filename,basename);
	strcat(filename,".tr");
	FILE * ft = fopen(filename,"r");
	for(fila=0;fila<trep->part;fila++){
		for(columna=0;columna<trep->part;columna++){  

			rep=trep->submatrices[fila][columna];
			if(rep->numberOfEdges==0)
				continue;
			load(rep->bt,ft);
			rep->bt_len = rep->bt->n;

		}
	}
	fclose(ft);

	strcpy(filename,basename);
	strcat(filename,".lv");
	FILE * fn = fopen(filename,"r");
	for(fila=0;fila<trep->part;fila++){
		for(columna=0;columna<trep->part;columna++){  

			rep=trep->submatrices[fila][columna];

			if(rep->numberOfEdges==0)
				continue;
			load(rep->bn,fn);
			rep->bn_len = rep->bn->n;
		}
	}
	fclose(fn);  

	free(filename);
	return trep;
}


uint * compactTreeAdjacencyList(TREP * trep, int x){
	trep->info[0]=0;
	MREP * rep;
	uint j,i;

	uint xrelatIn = x/trep->tamSubm;//*trep->part;
	uint * listady;
	x = x%trep->tamSubm;
	

	for(i=0;i<trep->part;i++){
		rep=trep->submatrices[xrelatIn][i];

		trep->columna=i;
		listady = compactAdjacencyList(trep,rep,x);

	}
	return trep->info;

}





uint * compactTreeInverseList(TREP * trep, int y){
	trep->info[0]=0;
	MREP * rep;
	uint j,i;

	uint yrelatIn = y/trep->tamSubm;//*trep->part;
	uint * listady;
	y = y%trep->tamSubm;
	
	for(i=0;i<trep->part;i++){
		rep=trep->submatrices[i][yrelatIn];
		trep->fila= i;
		listady = compactInverseList(trep,rep,y);
	}
	return trep->info;

}


void partialdestroyTreeRepresentation(TREP *trep){
	int i,j;
	for(i=0;i<trep->part;i++){
		for(j=0;j<trep->part;j++)
			free(trep->submatrices[i][j]);
		free(trep->submatrices[i]);
	}
	free(trep->submatrices);
	free(trep->div_level_table1);
	free(trep->div_level_table2);

	free(trep->info2[0]);
	free(trep->info2[1]);
	free(trep->info);
	free(trep->element);
	free(trep->basex);
	free(trep->basey);

	free(trep);		
}




void destroyTreeRepresentation(TREP *trep){
	int i,j;
	for(i=0;i<trep->part;i++){
		for(j=0;j<trep->part;j++)
			destroyRepresentation(trep->submatrices[i][j]);
		free(trep->submatrices[i]);
	}
	free(trep->submatrices);
	free(trep->div_level_table1);
	free(trep->div_level_table2);

	free(trep->info2[0]);
	free(trep->info2[1]);
	free(trep->info);
	free(trep->element);
	free(trep->basex);
	free(trep->basey);

	free(trep->basep1);
	free(trep->basep2);
	free(trep->baseq1);
	free(trep->baseq2);

	free(trep->words);
	free(trep);		
}

void   compressInformationLeaves(TREP * trep){
	//Traverse the leaves, compute frequency of submatrices and encoded them with DACs
	MREP * rep;

	uint fila,columna;

	ulong totalLeaves=0;
	uint i,j;

	for(fila=0;fila<trep->part;fila++){
		for(columna=0;columna<trep->part;columna++){
			rep=trep->submatrices[fila][columna];
			if(trep->submatrices[fila][columna]->numberOfEdges!=0)
				totalLeaves +=rep->nleaves;	
		}
	}
	
	unsigned char * ilchar = (unsigned char *) malloc(sizeof(unsigned char)*totalLeaves*K2_3_char);
	
	unsigned char *aWord;
	unsigned int size;
	initialize(totalLeaves);
	//Creating the vocabulary
	uint ilpos=0;
	uint jj;
	for(fila=0;fila<trep->part;fila++){
		for(columna=0;columna<trep->part;columna++){
			rep=trep->submatrices[fila][columna];
			if(trep->submatrices[fila][columna]->numberOfEdges==0)
				continue;
			for(i=0;i<rep->nleaves;i++){
				
				//Using hash for the vocabulary...
				aWord=&(ilchar[ilpos]);  //the word parsed.
				
				for(j=0;j<K2_3;j++){
					if(bitget(rep->leavesInf,i*K2_3+j))
						bitsetchar(aWord,j);
					else
						bitcleanchar(aWord,j);
					
				}
				
				size = K2_3_char;
				j = search ((unsigned char *)aWord, size, &addrInTH );

				if (j==zeroNode) {
					insertElement ((unsigned char *) aWord, size, &addrInTH);
					hash[addrInTH].weight = 0;
					hash[addrInTH].size = 0;
					hash[addrInTH].len = K2_3_char;
					positionInTH[zeroNode] = addrInTH;
					zeroNode++;
				}

				hash[addrInTH].weight +=1;

				ilpos+=K2_3_char;		
			}
		}
	}		
	
	
	trep->zeroNode = zeroNode;
	trep->lenWords = K2_3_char;
	
	//Compressing leaves
	
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

	//Aplies qsort to the words with frequency > 1.
	qsort(positionInTH,k,sizeof(unsigned int),comparaFrecListaDesc);

	//Generates codes sequentially

	ulong totalLeavesCount=0;
	for(i=0;i<zeroNode;i++){
		hash[positionInTH[i]].codeword = i;
		totalLeavesCount += hash[positionInTH[i]].weight;
	}
	

	/********************** Beginning of the second pass **********************/
	
	uint icont=0;

	unsigned int codeword;
	unsigned int tam;
	

	//Compressing leaves submatrices...

	ilpos=0;
	for(fila=0;fila<trep->part;fila++){
		for(columna=0;columna<trep->part;columna++){
			rep=trep->submatrices[fila][columna];
			if(trep->submatrices[fila][columna]->numberOfEdges==0)
				continue;
			unsigned int tamTotal = 0;
			uint * listIL = (uint *) malloc(sizeof(uint)*rep->nleaves);
			uint listILCount =0;


			for(i=0;i<rep->nleaves;i++){

						aWord=&(ilchar[ilpos]);  //the word parsed.
						size = K2_3_char;
						j = search ((unsigned char *)aWord, size, &addrInTH );
						
						listIL[listILCount++]=hash[addrInTH].codeword;

						ilpos+=K2_3_char;		

					}



					rep->compressIL = createFT(listIL,rep->nleaves);

					free(listIL);
				}
			}
			trep->words = (unsigned char *) malloc(sizeof(unsigned char)*trep->zeroNode*trep->lenWords);

			int wc = 0;
			for (i=0;i<zeroNode;i++){
				for(j=0;j<trep->lenWords;j++){
					trep->words[wc++]=hash[positionInTH[i]].word[j];
				}
			}

			free(ilchar);

		}


		void  openPartialFile(TREP * trep, char * basename){

			char *filename = (char *) malloc(sizeof(char)*(strlen(basename)+5));
			strcpy(filename,basename);
			strcat(filename,".tr");
			_ftr = fopen(filename,"w");


			strcpy(filename,basename);
			strcat(filename,".lv");
			_flv = fopen(filename,"w");


			strcpy(filename,basename);
			strcat(filename,".il");
			_fil = fopen(filename,"w");

			fwrite(&(trep->part),sizeof(uint),1,_fil);
			fwrite(&(trep->tamSubm),sizeof(uint),1,_fil);

			fwrite(&(trep->numberOfNodes),sizeof(uint),1,_fil);
			fwrite(&(trep->numberOfEdges),sizeof(ulong),1,_fil);

			fwrite(&(trep->repK1),sizeof(uint),1,_fil);
			fwrite(&(trep->repK2),sizeof(uint),1,_fil);
			fwrite(&(trep->maxRealLevel1),sizeof(uint),1,_fil);
			fwrite(&(trep->maxLevel1),sizeof(uint),1,_fil);
			fwrite(&(trep->maxLevel2),sizeof(uint),1,_fil);

			free(filename);

		}		  


		void  closePartialFile(){
			fclose(_ftr);
			fclose(_flv);  

			fwrite(&(_ttlLeaves),sizeof(ulong),1,_fil);

			fclose(_fil);
		}		  


		void   partialSave(TREP * trep, MREP * rep, uint fila, uint columna){

			if(rep->numberOfEdges!=0)
				_ttlLeaves +=rep->nleaves;	


			if(rep->numberOfEdges!=0)
				save(rep->bt,_ftr);


			if(rep->numberOfEdges!=0)
				save(rep->bn,_flv);
			


			destroyBitRankW32Int(rep->bt);
			destroyBitRankW32Int(rep->bn);




	//Save IL without compression



			
			fwrite(&(rep->numberOfNodes),sizeof(uint),1,_fil);
			fwrite(&(rep->numberOfEdges),sizeof(ulong),1,_fil);
			if(rep->numberOfEdges==0)
				return;
			fwrite(&(rep->cutBt),sizeof(uint),1,_fil);
			fwrite(&(rep->lastBt1_len),sizeof(uint),1,_fil);
			fwrite(&(rep->nleaves),sizeof(uint),1,_fil);
			if(rep->numberOfEdges==0)
				return;
			fwrite(rep->leavesInf,sizeof(uint),(rep->nleaves*K2_2*K2_2+W-1)/W,_fil);
			free(rep->leavesInf);

		}





		void   saveBeforeCompressInformationLeaves(TREP * trep, char * basename){
	//Recorrer todas las ils, calcular frecuencias y usar FT*
			MREP * rep;

			uint fila,columna;

			ulong totalLeaves=0;
			uint i,j;

			for(fila=0;fila<trep->part;fila++){
				for(columna=0;columna<trep->part;columna++){
					rep=trep->submatrices[fila][columna];
					if(trep->submatrices[fila][columna]->numberOfEdges!=0)
						totalLeaves +=rep->nleaves;	
				}
			}



	/*PRIMERA PARTE DE SAVETREE*/
			char *filename = (char *) malloc(sizeof(char)*(strlen(basename)+5));
			strcpy(filename,basename);
			strcat(filename,".tr");
			FILE * ft = fopen(filename,"w");
			uint part= trep->part;


			for(fila=0;fila<part;fila++){
				for(columna=0;columna<part;columna++){
					if(trep->submatrices[fila][columna]->numberOfEdges!=0)
						save(trep->submatrices[fila][columna]->bt,ft);
				}
			}
			fclose(ft);


			strcpy(filename,basename);
			strcat(filename,".lv");
			FILE * fn = fopen(filename,"w");
			for(fila=0;fila<part;fila++){
				for(columna=0;columna<part;columna++){
					if(trep->submatrices[fila][columna]->numberOfEdges!=0)
						save(trep->submatrices[fila][columna]->bn,fn);
				}
			}
			fclose(fn);  









		/*Primera parte de destroyTreeRepresentation*/

			for(i=0;i<trep->part;i++){
				for(j=0;j<trep->part;j++){
					destroyBitRankW32Int(trep->submatrices[i][j]->bt);
					destroyBitRankW32Int(trep->submatrices[i][j]->bn);}
				}
				free(trep->div_level_table1);
				free(trep->div_level_table2);

				free(trep->info2[0]);
				free(trep->info2[1]);
				free(trep->info);
				free(trep->element);
				free(trep->basex);
				free(trep->basey);


		/*FIn de destroy*/



	//GRABAR IL sin comprimir


				strcpy(filename,basename);
				strcat(filename,".il");
				FILE * fvil = fopen(filename,"w");

				fwrite(&(trep->part),sizeof(uint),1,fvil);
				fwrite(&(trep->tamSubm),sizeof(uint),1,fvil);

				fwrite(&(trep->numberOfNodes),sizeof(uint),1,fvil);
				fwrite(&(trep->numberOfEdges),sizeof(ulong),1,fvil);

				fwrite(&(trep->repK1),sizeof(uint),1,fvil);
				fwrite(&(trep->repK2),sizeof(uint),1,fvil);
				fwrite(&(trep->maxRealLevel1),sizeof(uint),1,fvil);
				fwrite(&(trep->maxLevel1),sizeof(uint),1,fvil);
				fwrite(&(trep->maxLevel2),sizeof(uint),1,fvil);


				fwrite(&(totalLeaves),sizeof(ulong),1,fvil);

				for(fila=0;fila<part;fila++){
					for(columna=0;columna<part;columna++){  
						rep=trep->submatrices[fila][columna];



						fwrite(&(rep->numberOfNodes),sizeof(uint),1,fvil);
						fwrite(&(rep->numberOfEdges),sizeof(ulong),1,fvil);
						if(trep->submatrices[fila][columna]->numberOfEdges==0)
							continue;
						fwrite(&(rep->cutBt),sizeof(uint),1,fvil);
						fwrite(&(rep->lastBt1_len),sizeof(uint),1,fvil);

						fwrite(&(rep->nleaves),sizeof(uint),1,fvil);

						if(trep->submatrices[fila][columna]->numberOfEdges==0)
							continue;
						fwrite(rep->leavesInf,sizeof(uint),(rep->nleaves*K2_2*K2_2+W-1)/W,fvil);

					}
				}

				fwrite(&(totalLeaves),sizeof(ulong),1,fvil);


				unsigned char *aWord;
				unsigned int size;
				unsigned char * buffer;

				buffer = (unsigned char*)malloc(sizeof(unsigned char)*K2_3_char);
				ulong contadorTL = 0;
				for(fila=0;fila<trep->part;fila++){
					for(columna=0;columna<trep->part;columna++){
						rep=trep->submatrices[fila][columna];
						if(trep->submatrices[fila][columna]->numberOfEdges==0)
							continue;
						for(i=0;i<rep->nleaves;i++){
				//Usando una hash para el vocabulario...
				aWord=buffer;  //the word parsed.
				
				for(j=0;j<K2_3;j++){
					if(bitget(rep->leavesInf,i*K2_3+j))
						bitsetchar(aWord,j);
					else
						bitcleanchar(aWord,j);
					
				}
				size= K2_3_char;
				//Guardando
				fwrite(aWord,sizeof(char),size,fvil);
				contadorTL++;
				
			}
		}
	}
	fclose(fvil);   
				//FIN GUARDADO

	
	
	unsigned char * ilchar = (unsigned char *) malloc(sizeof(unsigned char)*totalLeaves*K2_3_char);
	
	
	initialize(totalLeaves);
	//Creación del vocabulario
	uint ilpos=0;
	uint jj;
	for(fila=0;fila<trep->part;fila++){
		for(columna=0;columna<trep->part;columna++){
			rep=trep->submatrices[fila][columna];
			if(trep->submatrices[fila][columna]->numberOfEdges==0)
				continue;
			for(i=0;i<rep->nleaves;i++){
				
				//Usando una hash para el vocabulario...
				aWord=&(ilchar[ilpos]);  //the word parsed.
				
				for(j=0;j<K2_3;j++){
					if(bitget(rep->leavesInf,i*K2_3+j))
						bitsetchar(aWord,j);
					else
						bitcleanchar(aWord,j);
					
				}

				size= K2_3_char;
				j = search ((unsigned char *)aWord, size, &addrInTH );

				if (j==zeroNode) {
					insertElement ((unsigned char *) aWord, size, &addrInTH);
					hash[addrInTH].weight = 0;
					hash[addrInTH].size = 0;
					hash[addrInTH].len = K2_3_char;
					positionInTH[zeroNode] = addrInTH;
					zeroNode++;
				}

				hash[addrInTH].weight +=1;

				ilpos+=K2_3_char;		
			}
		}
	}		
	
	
	trep->zeroNode = zeroNode;
	trep->lenWords = K2_3_char;
	
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

	//Aplies qsort to the words with frequency > 1.
	qsort(positionInTH,k,sizeof(unsigned int),comparaFrecListaDesc);

	//Generates codes sequentially

	/* Compresion ETDC
	GeneraCodigosETDC (zeroNode);
	*/
	
	ulong totalLeavesCount=0;
	for(i=0;i<zeroNode;i++){
		hash[positionInTH[i]].codeword = i;
		totalLeavesCount += hash[positionInTH[i]].weight;
	}
	

	/********************** Beginning of the second pass **********************/
	
	uint icont=0;

	unsigned int codeword;
	unsigned int tam;
	

				//Compactando información de las hojas...

	ilpos=0;
	for(fila=0;fila<trep->part;fila++){
		for(columna=0;columna<trep->part;columna++){
			rep=trep->submatrices[fila][columna];
			if(trep->submatrices[fila][columna]->numberOfEdges==0)
				continue;
			unsigned int tamTotal = 0;
			uint * listIL = (uint *) malloc(sizeof(uint)*rep->nleaves);
			uint listILCount =0;


			for(i=0;i<rep->nleaves;i++){
				aWord=&(ilchar[ilpos]);  //the word parsed.
				size = K2_3_char;
				j = search ((unsigned char *)aWord, size, &addrInTH );

				listIL[listILCount++]=hash[addrInTH].codeword;

				ilpos+=K2_3_char;		

			}



			rep->compressIL = createFT(listIL,rep->nleaves);

			free(listIL);
		}
	}
	free(ilchar);

	trep->words = (unsigned char *) malloc(sizeof(unsigned char)*trep->zeroNode*trep->lenWords);

	int wc = 0;
	for (i=0;i<zeroNode;i++){
		for(j=0;j<trep->lenWords;j++){
			trep->words[wc++]=hash[positionInTH[i]].word[j];
		}
	}



	/* SEGUNDA PARTE DE SAVETREE*/

	strcpy(filename,basename);
	strcat(filename,".voc");
	FILE * fv = fopen(filename,"w");

	fwrite(&(trep->part),sizeof(uint),1,fv);
	fwrite(&(trep->tamSubm),sizeof(uint),1,fv);

	fwrite(&(trep->numberOfNodes),sizeof(uint),1,fv);
	fwrite(&(trep->numberOfEdges),sizeof(ulong),1,fv);

	fwrite(&(trep->repK1),sizeof(uint),1,fv);
	fwrite(&(trep->repK2),sizeof(uint),1,fv);
	fwrite(&(trep->maxRealLevel1),sizeof(uint),1,fv);
	fwrite(&(trep->maxLevel1),sizeof(uint),1,fv);
	fwrite(&(trep->maxLevel2),sizeof(uint),1,fv);

		//  fwrite(&(trep->maxLevel),sizeof(uint),1,fv);
		//  
		  fwrite(&trep->zeroNode,sizeof(uint),1,fv); //stores the number of words of the vocabulary
		  fwrite(&trep->lenWords,sizeof(uint),1,fv);
			//Writes the vocabulary to disk.

		  for (i=0;i<zeroNode;i++)
		  	fwrite(hash[positionInTH[i]].word,sizeof(char),trep->lenWords,fv);

		  fclose(fv);
		  
		  strcpy(filename,basename);
		  strcat(filename,".cil");
		  FILE * fi = fopen(filename,"w");


		  for(fila=0;fila<part;fila++){
		  	for(columna=0;columna<part;columna++){  
		  		rep=trep->submatrices[fila][columna];

		  		fwrite(&(rep->numberOfNodes),sizeof(uint),1,fi);
		  		fwrite(&(rep->numberOfEdges),sizeof(ulong),1,fi);
		  		if(trep->submatrices[fila][columna]->numberOfEdges==0)
		  			continue;
		  		fwrite(&(rep->cutBt),sizeof(uint),1,fi);
		  		fwrite(&(rep->lastBt1_len),sizeof(uint),1,fi);
		  		fwrite(&(rep->nleaves),sizeof(uint),1,fi);
		  		saveFT(rep->compressIL,fi);
		  	}
		  }
		  fclose(fi);   





		  free(filename);





}



uint compactTreeCheckLink(TREP * trep, uint x, uint y){
	int i,ii,j,jj;
	uint k, divlevel, p,q, pnew,qnew;


	int K1K1 = K1*K1;
	int K2K2 = K2*K2;

	MREP * rep;

	i=x/trep->tamSubm;
	j=y/trep->tamSubm;
	uint factorAdjust = K1K1/K2K2;

	p = x % trep->tamSubm;
	q = y % trep->tamSubm;

	rep=trep->submatrices[i][j];

	if(rep->numberOfEdges==0)
		return 0;



	uint nleaf,posInf,realvalue,  nleafrelat;
	uint prelat, qrelat;
	int node,div_level,xrelat;

	prelat = p;
	qrelat = q;

	node = 0;

	for(i=0;i<trep->maxLevel1-1;i++){

		div_level = trep->div_level_table1[i];

		node +=prelat/div_level*K1 + qrelat/div_level;

		prelat = prelat%div_level;
		qrelat = qrelat%div_level;


		if(isBitSet(rep->bt,node)==0)
			return 0;

		node = 	rank(rep->bt,node)*K1K1;

	}


	for(i=trep->maxLevel1-1;i<trep->maxLevel1;i++){

		div_level = trep->div_level_table1[i];

		node +=prelat/div_level*K1 + qrelat/div_level;
		prelat = prelat%div_level;
		qrelat = qrelat%div_level;


		if(isBitSet(rep->bt,node)==0)
			return 0;

		node = rep->cutBt + (rank(rep->bt,node)*K1K1-rep->cutBt)/factorAdjust;
	}

	uint cutPreRank = rank(rep->bt,rep->lastBt1_len-1);



	for(i=0;i<trep->maxLevel2-1;i++){

		div_level = trep->div_level_table2[i];

		node+=prelat/div_level*K2 + qrelat/div_level;

		prelat = prelat%div_level;
		qrelat = qrelat%div_level;


		if(isBitSet(rep->bt,node)==0)
			return 0;

		node =rep->cutBt + (rank(rep->bt,node-1)-cutPreRank)*K2K2;

	}


	nleaf = node-rep->bt_len;

	nleafrelat = nleaf + (prelat/K2_2)*K2+qrelat/K2_2;




	if(isBitSet(rep->bn,nleafrelat)){

		posInf = rank(rep->bn,nleafrelat);

		realvalue = accessFT(rep->compressIL,posInf);				
		if(bitgetchar(&(trep->words[realvalue*trep->lenWords]),(qrelat%K2_2+(prelat%K2_2)*K2_2)))
			return 1;
	}

	return 0;	

}




void compactRangeQuery(TREP * trep, MREP * rep,uint p1, uint p2, uint q1, uint q2){
	int i,ii,j;
	uint y,k, divlevel, p1new, p2new, q1new, q2new;
	uint realvalue, leaf, posInf;


	int K1K1 = K1*K1;
	int K2K2 = K2*K2;
	uint AJUSTEK1K2 = K1K1/K2K2;
	if(rep->numberOfEdges==0)
		return;    

	trep->iniq=-1;
	trep->finq=-1;
	uint nleaf, nleafrelat, dp, dq;

	uint totalAdyNodes =0;
	int  queuecont, conttmp,node, node2;

	AddItem3(trep,0,0,0,p1,p2,q1,q2);

	queuecont = 1;
	for(ii=0;ii<trep->maxLevel1-1;ii++){
		conttmp = 0;

		divlevel = trep->div_level_table1[ii];

		for(k=0;k<queuecont;k++){
			node = trep->element[trep->iniq];
			p1= trep->basep1[trep->iniq];
			p2= trep->basep2[trep->iniq];
			q1= trep->baseq1[trep->iniq];
			q2= trep->baseq2[trep->iniq];
			dp = trep->basex[trep->iniq];
			dq = trep->basey[trep->iniq];

			for(i=p1/divlevel;i<=p2/divlevel;i++){
				if(i==p1/divlevel)
					p1new=p1 % divlevel;
				else
					p1new=0;
				if(i==p2/divlevel)
					p2new=p2 % divlevel;
				else
					p2new=divlevel-1;

				for(j=q1/divlevel;j<=q2/divlevel;j++){
					if(j==q1/divlevel)
						q1new=q1 % divlevel;
					else
						q1new=0;
					if(j==q2/divlevel)
						q2new=q2 % divlevel;
					else
						q2new=divlevel-1;    

					node2= node+K1*i+j;
                    //fprintf(stderr,"node: %d, node2: %d\n",node, node2);
										//fprintf(stderr,"en bt_len: %d, accediendo a %d\n",rep->bt_len, node2);
					if(isBitSet(rep->bt,node2)){

						conttmp++;

						AddItem3(trep,rank(rep->bt,node2)*K1K1,dp + divlevel*i,dq + divlevel*j,p1new,p2new,q1new,q2new);
					}

				}
			}

			RemoveItem3(trep);
		}
		queuecont = conttmp;
	}




	for(ii=trep->maxLevel1-1;ii<trep->maxLevel1;ii++){
		conttmp = 0;


		divlevel = trep->div_level_table1[ii];

		for(k=0;k<queuecont;k++){
			node = trep->element[trep->iniq];
			p1= trep->basep1[trep->iniq];
			p2= trep->basep2[trep->iniq];
			q1= trep->baseq1[trep->iniq];
			q2= trep->baseq2[trep->iniq];
			dp = trep->basex[trep->iniq];
			dq = trep->basey[trep->iniq];


			for(i=p1/divlevel;i<=p2/divlevel;i++){
				if(i==p1/divlevel)
					p1new=p1 % divlevel;
				else
					p1new=0;
				if(i==p2/divlevel)
					p2new=p2 % divlevel;
				else
					p2new=divlevel-1;

				for(j=q1/divlevel;j<=q2/divlevel;j++){
					if(j==q1/divlevel)
						q1new=q1 % divlevel;
					else
						q1new=0;
					if(j==q2/divlevel)
						q2new=q2 % divlevel;
					else
						q2new=divlevel-1;    

					node2= node+K1*i+j;

					if(isBitSet(rep->bt,node2)){

						conttmp++;


						AddItem3(trep,rep->cutBt + (rank(rep->bt,node2)*K1K1-rep->cutBt)/AJUSTEK1K2,dp + divlevel*i,dq + divlevel*j,p1new,p2new,q1new,q2new);

					}
				}
			}

			RemoveItem3(trep);
		}
		queuecont = conttmp;
	}

	uint cutPreRank = rank(rep->bt,rep->lastBt1_len-1);



	for(ii=0;ii<trep->maxLevel2-1;ii++){
		conttmp = 0;


		divlevel = trep->div_level_table2[ii];

		for(k=0;k<queuecont;k++){
			node = trep->element[trep->iniq];
			p1= trep->basep1[trep->iniq];
			p2= trep->basep2[trep->iniq];
			q1= trep->baseq1[trep->iniq];
			q2= trep->baseq2[trep->iniq];
			dp = trep->basex[trep->iniq];
			dq = trep->basey[trep->iniq];


			for(i=p1/divlevel;i<=p2/divlevel;i++){
				if(i==p1/divlevel)
					p1new=p1 % divlevel;
				else
					p1new=0;
				if(i==p2/divlevel)
					p2new=p2 % divlevel;
				else
					p2new=divlevel-1;

				for(j=q1/divlevel;j<=q2/divlevel;j++){
					if(j==q1/divlevel)
						q1new=q1 % divlevel;
					else
						q1new=0;
					if(j==q2/divlevel)
						q2new=q2 % divlevel;
					else
						q2new=divlevel-1;    


					node2= node+K2*i+j;

					if(isBitSet(rep->bt,node2)){

						conttmp++;

						AddItem3(trep,rep->cutBt + (rank(rep->bt,node2-1)-cutPreRank)*K2K2,dp + divlevel*i,dq + divlevel*j,p1new,p2new,q1new,q2new);

					}
				}
			}

			RemoveItem3(trep);
		}
		queuecont = conttmp;
	}


	while(trep->iniq<=trep->finq){
		nleaf = trep->element[trep->iniq]-rep->bt_len;
		p1= trep->basep1[trep->iniq];
		p2= trep->basep2[trep->iniq];
		q1= trep->baseq1[trep->iniq];
		q2= trep->baseq2[trep->iniq];
		dp = trep->basex[trep->iniq];
		dq = trep->basey[trep->iniq];

		int jj;

		for(i=p1/K2_2;i<=p2/K2_2;i++){
			if(i==p1/K2_2)
				p1new=p1 % K2_2;
			else
				p1new=0;
			if(i==p2/K2_2)
				p2new=p2 % K2_2;
			else
				p2new=K2_2-1;

			for(j=q1/K2_2;j<=q2/K2_2;j++){

				if(j==q1/K2_2)
					q1new=q1 % K2_2;
				else
					q1new=0;
				if(j==q2/K2_2)
					q2new=q2 % K2_2;
				else
					q2new=K2_2-1;
				leaf = nleaf + i*K2+j;

				if(isBitSet(rep->bn,leaf)){


					posInf= rank(rep->bn,leaf);
					realvalue = accessFT(rep->compressIL,posInf);    


					for(ii=p1new;ii<=p2new;ii++)
						for(jj=q1new;jj<=q2new;jj++)

							if(bitgetchar(&(trep->words[realvalue*trep->lenWords]),(jj%K2_2+(ii%K2_2)*K2_2))){
								trep->info2[0][0]++;
								trep->info2[0][trep->info2[0][0]]=dp+i*K2_2+ii+trep->fila*trep->tamSubm;
								trep->info2[1][trep->info2[0][0]]=dq+j*K2_2+jj+trep->columna*trep->tamSubm;
							}
						}
					}
				}
				RemoveItem3(trep);
			}


			return;

		}

uint ** compactTreeRangeQuery(TREP * trep, uint p1, uint p2, uint q1, uint q2){

	trep->info2[0][0]=0;

	MREP * rep;
	uint j,i ,p1new, p2new, q1new, q2new;

	uint * listady;


	for(i=p1/trep->tamSubm;i<=p2/trep->tamSubm;i++){
		if(i==p1/trep->tamSubm)
			p1new=p1 % trep->tamSubm;
		else
			p1new=0;
		if(i==p2/trep->tamSubm)
			p2new=p2 % trep->tamSubm;
		else
			p2new=trep->tamSubm-1;
		for(j=q1/trep->tamSubm;j<=q2/trep->tamSubm;j++){
			if(j==q1/trep->tamSubm)
				q1new=q1 % trep->tamSubm;
			else
				q1new=0;
			if(j==q2/trep->tamSubm)
				q2new=q2 % trep->tamSubm;
			else
				q2new=trep->tamSubm-1;



			rep=trep->submatrices[i][j];

			trep->fila=i;
			trep->columna=j;

			compactRangeQuery(trep,rep, p1new,p2new,q1new,q2new);

		}
	}
	return trep->info2;

}

