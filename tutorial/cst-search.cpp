#include <sdsl/suffix_trees.hpp>
#include <iostream>
#include <string>

using namespace sdsl;
using namespace std;


template<class t_cst, class t_pat>
void execute(const char* input, uint8_t num_bytes, t_pat& pat, const char* format, char sentinel)
{
    typedef typename t_cst::node_type node_t;
    t_cst cst;
    construct_im(cst, input, num_bytes);
    csXprintf(cout, format, cst, sentinel);

    cout << "pattern \"" << pat << "\"" << endl;
    cout << "---- backward search step by step ----" << endl;
    {
        uint64_t lb=0, rb=cst.size()-1;
        for (auto it=pat.end(); it != pat.begin() and lb <= rb;) {
            --it;
            if (backward_search(cst.csa, lb, rb, (typename t_cst::char_type)*it, lb, rb) > 0) {
                cout << "[" << lb << "," << rb << "]" << endl;
                cout << "matched " << *it << endl;
            }
        }
    }

    cout << "---- backward search for the pattern  ----" << endl;
    {
        uint64_t lb=0, rb=cst.size()-1;
        backward_search(cst.csa, lb, rb, pat.begin(), pat.end(), lb, rb);
        cout << "size = " << rb+1-lb << endl;
    }
    cout << "---- count pattern occurrences  ----" << endl;
    {
        cout << "count(cst.csa, \"" << pat << "\")=" << count(cst.csa, pat.begin(), pat.end()) << endl;
    }
    cout << "---- locate the pattern  ----" << endl;
    {
        typedef int_vector<> vec_t;
        vec_t occ;
        cout << "locate(cst.csa, \"" << pat << "\")=" << locate(cst.csa, pat.begin(), pat.end(), occ) << endl;
        for (auto it=occ.begin(); it != occ.end(); ++it) {
            cout << *it << " ";
        }
        cout << endl;
    }
    cout << "---- extract text  ----" << endl;
    {
        cout << "extract(csa,0,csa.size())=\"" << extract(cst.csa, 0, cst.csa.size()-1) << "\"" << endl;
    }

    cout << "---- forward search step by step  ----" << endl;
    {
        node_t v = cst.root();
        auto it = pat.begin();
        for (uint64_t char_pos=0; it != pat.end(); ++it) {
            if (forward_search(cst, v, it-pat.begin(), *it, char_pos) > 0) {
                cout << it-pat.begin() << "-[" << cst.lb(v) << "," << cst.rb(v) << "]" << endl;
                cout << "matched " << *it << endl;
            } else {
                break;
            }
        }
    }
    cout << "---- count pattern occurrences ----" << endl;
    {
        cout << "count(cst, \"" << pat << "\")=" << count(cst, pat.begin(), pat.end()) << endl;
    }
    cout << "---- extract text  ----" << endl;
    {
        cout << "extract(cst,cst.select_leaf(cst.csa.isa[0]+1))=\""
             << extract(cst, cst.select_leaf(cst.csa.isa[0]+1)) << "\"" << endl;
    }

}


int main()
{
    {
        cout << "# Byte alphabet example\n" << endl;
        string pat("brac");
        execute<cst_sct3<> >("abracadabra#bracda", 1, pat, "%2I %3S %T", '$');
        cout << "\n\n" << endl;
    }

    {
        cout << "# Integer alphabet example\n" << endl;
        int_vector<> pat(2);
        pat[0] = 801; pat[1] = 444;
        execute<cst_sct3<csa_bitcompressed<int_alphabet<> > > >("2 801 543 293 597 801 444 444 293", 'd', pat, "%2I %3S %:4T", '0');
    }
}
