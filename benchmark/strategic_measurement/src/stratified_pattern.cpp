#include "types.hpp"
#include <cstdlib>
#include <iostream>
#include <queue>
#include <string>

using namespace sdsl;
using namespace std;

typedef pair<uint64_t,uint64_t> pii_t;
typedef priority_queue<pii_t> pq_t;

int main(int argc, char* argv[])
{
    if (argc < 8) {
        cout << "Usage: ./" << argv[0] << " file tmp_dir LEN MIN_OCC MAX_OCC RES_LEN pat_file" << endl;
        cout << " Calculates the average and medium number of occurrences" << endl;
        cout << " of all substrings of length m in `file`." << endl;
        cout << " Lengths m are read from stdin. For each m we output:" << endl;
        cout << " TC_ID   = `TC_ID`" << endl;
        cout << " m       = `m`" << endl;
        cout << " avg_occ = `avg_occ`" << endl;
        cout << " med_occ = `med_occ`" << endl;
        return 1;
    }
    std::string tmp_dir;
    uint64_t len     = atoll(argv[3]);
    uint64_t min_occ = atoll(argv[4]);
    uint64_t max_occ = atoll(argv[5]);
    uint64_t res_len = atoll(argv[6]);
    cst_t cst;
    cache_config cconfig(false, argv[2], util::basename(argv[1]));
    if (!cache_file_exists("stratcst", cconfig)) {
        construct(cst, argv[1], cconfig, 1);
        store_to_cache(cst, "stratcst", cconfig);
    } else {
        load_from_cache(cst, "stratcst", cconfig);
    }

    pq_t pq; // max priority queue for potential pattern

    bit_vector short_leaf(cst.size(),0);
    for (uint64_t x = cst.size()-len+2; x < cst.size(); ++x) {
        short_leaf[cst.csa(x)] = 1;
    }
    uint64_t candidates = 0;
    // do a DFS-traversal through the tree
    for (cst_t::const_iterator it = cst.begin(); it != cst.end(); ++it) {
        if (it.visit() == 1) {  // first visit of the node
            cst_t::node_type v = *it;
            if (cst.size(v) < min_occ) {
                it.skip_subtree();
            } else { // cst.size(v) >= min_occ
                if ((cst.is_leaf(v) and !short_leaf[cst.lb(v)]) or cst.depth(v) >= len) {  //  >= len
                    if (cst.size(v) <= max_occ) {
                        ++candidates;
                        uint64_t p = cst.lb(v);
                        uint64_t r = lrand48()%cst.size();// generate rand value in [0..cst.size())
                        if (pq.size() < res_len) {
                            pq.push(pii_t(r, p));
                        } else if (pq.top().first > r) {  // now pq.size() == res_len
                            pq.pop(); // pop top most element
                            pq.push(pii_t(r, p));
                        }
                    }
                    it.skip_subtree();
                }
            }
        }
    }

    ofstream out(argv[7]);
    out<<"# number="<<pq.size()<<" length="<<len<<" file="<<argv[1];
    out<<" candidates="<<candidates<<std::endl;
    while (!pq.empty()) {
        uint64_t start_idx = cst.csa[ pq.top().second ];
        pq.pop();
        out << extract<string>(cst.csa, start_idx, start_idx+len-1);
    }
}
