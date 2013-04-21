#include <sdsl/suffix_trees.hpp>
#include <iostream>
#include <map>

using namespace sdsl;
using namespace std;

typedef cst_sada<csa_wt<>, lcp_dac<> > cst_t;
typedef map<uint64_t, uint64_t> muu_t;

int main(int argc, char* argv[])
{
    if (argc < 4) {
        cout << "Usage: ./" << argv[0] << " file tmp_dir TC_ID" << endl;
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
    cst_t cst;
    cache_config cconfig(false, argv[2], util::basename(argv[1]));
    construct(cst, argv[1], cconfig, 1);

    uint64_t m = 0;
    while (cin>> m) {   // we assume m >= cst.size()-1
        bit_vector short_leaf(cst.size(),0);
        for (uint64_t x = cst.size()-m+2; x < cst.size(); ++x) {
            short_leaf[cst.csa(x)] = 1;
        }
        uint64_t nocc = 0, npat = 0;
        muu_t occ_map;
        // do a DFS-traversal through the tree
        for (cst_t::const_iterator it = cst.begin(); it != cst.end(); ++it) {
            if (it.visit() == 1) {  // first visit of the node
                cst_t::node_type v = *it;
                if ((cst.is_leaf(v) and !short_leaf[cst.lb(v)]) or cst.depth(v) >= m) {    // if not >= m
                    nocc += cst.size(v)*cst.size(v);     // add occurrences
                    npat += cst.size(v);
                    occ_map[cst.size(v)]+=cst.size(v); // add to map for medium calculation
                    it.skip_subtree();      // skip subtree
                }
            }
        }
        uint64_t med = 0, pat = 0;
        for (muu_t::const_iterator it = occ_map.begin(); pat <= npat/2 and it != occ_map.end(); ++it) {
            pat += it->second;
            med = it->first;
        }
        cout << "# TC_ID   = " << argv[3] << endl;
        cout << "# m       = " << m << endl;
        cout << "# avg_occ = " << ((long double)nocc)/npat << endl;
        cout << "# med_occ = " << med << endl;
        cout << "npat ="<<npat<<endl;
    }
}
