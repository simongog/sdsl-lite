#include <stdio.h>
#include <math.h>
#include <string.h>
#include "kTree.h"



int main(int argc, char* argv[]){

	if(argc<2){
		fprintf(stderr,"USAGE: %s <GRAPH>\n",argv[0]);
		return(-1);
	}

	char *filename = (char *)malloc(sizeof(char)*(strlen(argv[1])+5));
	TREP * rep = loadTreeRepresentation(argv[1]);
	
  strcpy(filename,argv[1]);
  strcat(filename,".rrb");
	FILE *fr = fopen(filename,"w");
	fwrite(&(rep->numberOfNodes),sizeof(uint),1,fr);
  fwrite(&(rep->numberOfEdges),sizeof(ulong),1,fr);

	uint * listady, node, i, j,k, edge;
  //Using range queries
  	uint ** listady2;
	for(i=0;i<rep->numberOfNodes;i++){
  	node = -(i+1);
  	fwrite(&node,sizeof(uint),1,fr);
		uint ** resp = compactTreeRangeQuery(rep, i,i,0,rep->numberOfNodes);
		for(k=0;k<resp[0][0];k++)	{
		  	edge = resp[1][k+1] + 1;
		  	fwrite(&edge,sizeof(uint),1,fr);

  	}
  }
  
  
  fclose(fr);
  
  destroyTreeRepresentation(rep);
  free(filename);
  
  return 0;
}


