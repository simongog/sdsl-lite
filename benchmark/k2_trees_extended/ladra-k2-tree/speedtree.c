#include <stdio.h>
#include <math.h>
#include <string.h>
#include "kTree.h"
#include <stdbool.h>

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
		fprintf(stderr,"USAGE: %s <GRAPH> <QUERIES> (comp?)\n",argv[0]);
		return(-1);
	}

	bool comp = false;
	if (argc>2)
		comp = true;
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
    uint *l  = compactTreeAdjacencyList(trep, qry[i]);
    recovered += l[0];
  }
  t += stop_clock(); 
  t *= 1000; // to milliseconds

	fclose(list_fp);
	if (!comp){
		fprintf(stderr,"# adj_time = 0\n");
		fprintf(stderr,"# adj_check = 0\n");
		fprintf(stderr,"# neighbors_time = %f\n",t/recovered*1000);
		fprintf(stderr,"# neighbors_check = %lld\n",recovered);
	} else {
		fprintf(stderr,"# adj_time_comp = 0\n");
		fprintf(stderr,"# adj_check_comp = 0\n");
		fprintf(stderr,"# neighbors_time_comp = %f\n",t/recovered*1000);
		fprintf(stderr,"# neighbors_check_comp = %lld\n",recovered);
	}
	
 // destroyTreeRepresentation(trep);

  
  return 0;
}


