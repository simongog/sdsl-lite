#include <stdio.h>
#include <math.h>
#include <string.h>
#include "kTree.h"
#include <sys/time.h>
#include <sys/resource.h>


//static struct rusage resource_usage;

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
	FILE *f;
	uint nodes; 
	ulong edges;
	register ulong i;
	TREP * trep;

	if(argc<6){
		fprintf(stderr,"USAGE: %s <GRAPH> <name> <K1> <K2> [<max level K1> <S> <fastConstruction>]\n",argv[0]);
		return(-1);
	}
	
	fprintf(stderr, "Opening: %s\n", argv[1]);

	f = fopen(argv[1],"r");
	fread(&nodes,sizeof(uint),1,f);

	int _K1 = atoi(argv[3]);
	int _K2 = atoi(argv[4]);

	uint s=20;
	if(argc>6)
		s =  (atoi(argv[6]));
	else 
		s =  (atoi(argv[5]));
	uint tamSubm = 1<<s;


	if(argc>7){
		double t = 0;
		ticks= (double)sysconf(_SC_CLK_TCK);
  		start_clock();
		trep = compactCreateKTree(argv[1],argv[2],tamSubm, _K1,_K2, atoi(argv[5]));
		t += stop_clock(); 
  		t *= 1000; // to milliseconds
  		fprintf(stderr,"# constructs_time = %f\n",t/1000);

		struct rusage resource_usage;
    	getrusage(RUSAGE_SELF, &resource_usage);  		

  		fprintf(stderr,"# constructs_space = %ld\n",resource_usage.ru_maxrss*1024);
  		fprintf(stderr,"# constructs_space_vmem = 0\n");
  		fprintf(stderr,"# construct_morton_duration = 0\n");
  		fprintf(stderr,"# construct_bv_complete_duration = 0\n");
  		fprintf(stderr,"# construct_sort_duration = 0\n");
  		fprintf(stderr,"# construct_duration = 0\n");
  		fprintf(stderr,"# buildvec_duration = 0\n");
  		fprintf(stderr,"# subtree_constructor_duration = 0\n");
  		fprintf(stderr,"# constructor_call_duration = 0\n");
	}
	else{
		uint part;
		part = nodes/tamSubm+1;


		uint nodesOrig = nodes;
		nodes= tamSubm; 




		int max_level1;
		int max_real_level1;
		max_real_level1 = ceil(log(nodes)/log(_K1))-1;
		if(argc<7)
			max_level1 = ceil(log(nodes)/log(_K1))-L;
		else
			max_level1 = atoi(argv[5]);


		int nodes2 = 1;
		for(i=0;i<max_real_level1+1;i++){
			nodes2 *=_K1;
		}



		for(i=0;i<max_level1;i++)
			nodes2 = ceil((double)nodes2/_K1);

		int	max_level2 = ceil(log(nodes2)/log(_K2))-1;  

		fread(&edges,sizeof(ulong),1,f);


		uint nodes_read=0;

		uint foo;
		ulong foo2;
		int id1, id2;
		NODE * tree;
		MREP * rep;

		int fila,columna;




		trep = createTreeRep(nodesOrig,edges,part,tamSubm, max_real_level1, max_level1, max_level2,_K1,_K2);
		fclose(f);

		openPartialFile(trep,argv[2]);

		ulong edges_read=0;
		for(fila=0;fila<part;fila++){

			for(columna=0;columna<part;columna++){
				ulong edges_sub=0;

				numberNodes=0;
				numberLeaves=0;
				numberTotalLeaves=0;


				f = fopen(argv[1],"r");

				fread(&foo,sizeof(uint),1,f);
				fread(&foo2,sizeof(ulong),1,f);

				tree = createKTree(_K1,_K2,max_real_level1,max_level1,max_level2);

				nodes_read=0;
				edges_read=0;

				for(i=0;i<nodesOrig+edges;i++) {
					int k;
					fread(&k,sizeof(uint),1,f);
					if(k<0) {
						nodes_read++;
					}
					else {
						k--;

						id1=nodes_read-1;
						id2=k;

						if((id1>=fila*tamSubm)&&(id1<(fila+1)*tamSubm)&&(id2>=columna*tamSubm)&&(id2<(columna+1)*tamSubm)){

							insertNode(tree,id1-fila*tamSubm,id2-columna*tamSubm);
							edges_sub++;
						}
						edges_read++;
					}
				}


				fclose(f);

				MREP * rep;
				rep = createRepresentation(tree,nodes,edges_sub);


				insertIntoTreeRep(trep, rep, fila, columna);

				partialSave(trep,rep,fila,columna);	
			}
		}
	}
	partialdestroyTreeRepresentation(trep);
	closePartialFile();

	return 0;
}


