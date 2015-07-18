#include <iostream>
#include <sdsl/lcp.hpp>
#include <sdsl/suffix_arrays.hpp>

using namespace std;
using namespace sdsl;

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " file" << endl;
        return 1;
    }
    string file = argv[1];

    cache_config cc(false); // do not delete temp files after csa construction
    csa_wt<> csa;
    construct(csa, file, 1);

    cc.delete_files = true; // delete temp files after lcp construction
    lcp_wt<> lcp;
    construct(lcp, file, 1);

    if (csa.size() < 1000) {
        cout << csa << endl;
        cout << "-------" << endl;
        cout << lcp << endl;
    }
}
