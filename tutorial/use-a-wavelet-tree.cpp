#include <sdsl/wavelet_trees.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "Usage: " << argv[1] << " file" << endl;
        return 1;
    }
    wt_huff<rrr_vector<63>> wt;
    construct(wt, argv[1], 1);

    cout << "wt.size()="<< wt.size() << endl;
    cout << "wt.sigma ="<< wt.sigma << endl;
    if (wt.size() > 0) {
        // access an element
        cout << "wt[0]=" << wt[0] << endl;
        // rank an element (exclude)
        uint64_t r = wt.rank(wt.size(), wt[0]);
        cout << "wt.rank(wt.size(), wt[0])=" << r  << endl;
        // select element ()
        cout << "wt.select(r, wt[0]) = " << wt.select(r, wt[0]) << endl;
    }
}
