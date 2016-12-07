#include <stdio.h>
#include <math.h>
#include <string.h>
#include "kTree.h"
#include <stdbool.h>

struct pair {
    uint first;
    uint second;
};

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
  	struct pair *qry = (struct pair *) malloc(sizeof(struct pair)*queries);
  	fread(qry,sizeof(struct pair),queries,list_fp);
  	fprintf(stderr,"Processing %d queries\n",queries);
  	ticks= (double)sysconf(_SC_CLK_TCK);
  	start_clock();
  	uint i,l;
  	for(i=0;i<queries;i++) {
	    l  = compactTreeCheckLink(trep, qry[i].first, qry[i].second);
    	recovered += l;
  	}
  	t += stop_clock(); 
  	t *= 1000; // to milliseconds

	fclose(list_fp);
	if (!comp){
		fprintf(stderr,"# adj_time = %f\n",t/recovered*1000);
		fprintf(stderr,"# adj_check = %lld\n",recovered);
	} else {
		fprintf(stderr,"# adj_time_comp = %f\n",t/recovered*1000);
		fprintf(stderr,"# adj_check_comp = %lld\n",recovered);
	}
	
 // destroyTreeRepresentation(trep);

  
  return 0;
}


