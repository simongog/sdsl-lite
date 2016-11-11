#include <stdio.h>
#include <math.h>
#include <string.h>
#include "kTree.h"



int main(int argc, char* argv[]){

	if(argc<7){
		fprintf(stderr,"USAGE: %s <GRAPH> <p> <p1> <p2> <q1> <q2>\n",argv[0]);
		return(-1);
	}

	
	TREP * trep = loadTreeRepresentation(argv[1]);
	fprintf(stderr,"Adjacency list of node %d:\n",atoi(argv[2]));
	uint i;
	
	
  uint * listady = compactTreeAdjacencyList(trep, atoi(argv[2]));

  
  fprintf(stderr,"Number of neighbours: %d\n",listady[0]);
  for(i=0;i<listady[0];i++)
  	fprintf(stderr,"%d\n",listady[i+1]);
	fprintf(stderr,"\n");

  fprintf(stderr,"------ INVERSES --------\n");


  listady = compactTreeInverseList(trep, atoi(argv[2]));
  
  fprintf(stderr,"Number of neighbours: %d\n",listady[0]);
  for(i=0;i<listady[0];i++)
  	fprintf(stderr,"%d\n",listady[i+1]);
  	fprintf(stderr,"\n");



  
  uint p1=atoi(argv[3]);
  uint p2=atoi(argv[4]);
  uint q1=atoi(argv[5]);
  uint q2=atoi(argv[6]);
  uint ** respuesta = compactTreeRangeQuery(trep, p1, p2, q1, q2);

  fprintf(stderr,"--------------\n");
  fprintf(stderr,"Range: [%d,%d]-[%d,%d], total of links %d\n",p1,p2,q1,q2,respuesta[0][0]);
 for(i=0;i<respuesta[0][0];i++)
        fprintf(stderr,"(%d,%d)\t",respuesta[0][i+1],respuesta[1][i+1]);
        fprintf(stderr,"\n");

  
  destroyTreeRepresentation(trep);

  
  return 0;
}


