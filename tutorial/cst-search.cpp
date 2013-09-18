#include <sdsl/suffix_trees.hpp>
#include <iostream>
#include <string>

using namespace sdsl;
using namespace std;


template<class t_cst, class t_pat=typename t_cst::string_type>
void execute(const char* input, uint8_t num_bytes, t_pat pat, const char* format)
{
    typedef typename t_cst::node_type node_t;
    t_cst cst;
    construct_im(cst, input, num_bytes);
    csXprintf(cout, format, cst);

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
        cout << "count(cst.csa, \"" << pat << "\")=" << count(cst.csa, pat) << endl;
    }
    cout << "---- locate the pattern  ----" << endl;
    {
        auto occs = locate(cst.csa, pat);
        cout << "locate(cst.csa, \"" << pat << "\")=" << occs.size() << endl;
        cout << occs << endl;
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
        cout << "count(cst, \"" << pat << "\")=" << count(cst, pat) << endl;
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
        execute<cst_sct3<> >("abracadabra#bracda", 1, string("brac"), "%2I %3S %T");
        cout << "\n\n" << endl;
    }

    {
        cout << "# Integer alphabet example\n" << endl;
        execute<cst_sct3<csa_bitcompressed<int_alphabet<> > > >("2 801 543 293 597 801 444 444 293", 'd', {801, 444}, "%2I %3S %:4T");
    }
}
