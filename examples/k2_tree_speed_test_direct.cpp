#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <complex>
#include <sdsl/k2_tree.hpp>
#include <sdsl/k2_tree_algorithm.hpp>
#include <sdsl/bit_vectors.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <sys/times.h>

using namespace sdsl;

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
        fprintf(stderr,"USAGE: %s <path to serialized k tree (set templates correctly!)> query_file (binary format)\n",argv[0]);
        return(-1);
    }

    //char *filename = (char *)malloc(sizeof(char)*20);

//    const uint8_t k = 4;
    //typedef k2_tree_hybrid<k,k,k,k, bit_vector, bit_vector,true> k2_rrr;
//    typedef k2_tree<k, bit_vector, bit_vector, true, 4> k2_rrr;
        const uint8_t k = 4;
        //typedef k2_tree_hybrid<4,5,2,8, bit_vector, bit_vector> k2_rrr;
        //typedef k2_tree_partitioned<4, k2_rrr, true> k2_part;

    typedef k2_tree_hybrid<4,5,2,8, bit_vector, bit_vector, false> k2_rrr;
    typedef k2_tree_partitioned<4, k2_rrr, true> k2_part;

    k2_part k2tree;
    std::string fileName = argv[1];
    load_from_file(k2tree, fileName);

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
    std::vector<uint32_t> result;
    for(i=0;i<queries;i++) {
            k2tree.direct_links2(qry[i], result);
            recovered+= result.size();
    }
    t += stop_clock();
    t *= 1000; // to milliseconds

    fclose(list_fp);
    fprintf(stderr,"Recovered Nodes: %lld\n",recovered);
    fprintf(stderr,"Queries: %d\n",queries);
    fprintf(stderr,"Total time(ms): %f",t);
    fprintf(stderr,"Time per query: %f\n",t/queries);
    fprintf(stderr,"Time per link: %f\n",t/recovered);

    // destroyTreeRepresentation(trep);


    return 0;
}
