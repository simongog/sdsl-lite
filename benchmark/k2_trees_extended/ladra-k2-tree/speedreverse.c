#include <stdio.h>
#include <math.h>
#include <string.h>
#include "kTree.h"

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
/* end Time meassuring */

int main(int argc, char* argv[]){

	if(argc<2){
		fprintf(stderr,"USAGE: %s <GRAPH>\n",argv[0]);
		return(-1);
	}

	//char *filename = (char *)malloc(sizeof(char)*20);
	

		TREP * trep = loadTreeRepresentation(argv[1]);
	
		char * list_file = argv[2];
		
		
		
		FILE * list_fp = fopen(list_file,"r");
	uint queries;
	fread(&queries, sizeof(uint), 1, list_fp);
	ulong recovered = 0;
	double t = 0;
  uint *qry = (uint *) malloc(sizeof(uint)*queries);
  fread(qry,sizeof(uint),queries,list_fp);
  fprintf(stderr,"Processing %d queries\n",queries);
  ticks= (double)sysconf(_SC_CLK_TCK);
  start_clock();
  uint i;
  for(i=0;i<queries;i++) {
    uint *l  = compactTreeInverseList(trep, qry[i]);
    recovered += l[0];
  }
  t += stop_clock(); 
  t *= 1000; // to milliseconds

	fclose(list_fp);
	fprintf(stderr,"Recovered Nodes: %d\n",recovered);
	fprintf(stderr,"Queries: %d\n",queries);
	fprintf(stderr,"Total time(ms): %f",t);
	fprintf(stderr,"Time per query: %f\n",t/queries);
	fprintf(stderr,"Time per link: %f\n",t/recovered);

  destroyTreeRepresentation(trep);

  
  return 0;
}


