#include <sdsl/suffix_arrays.hpp>
#include <iostream>

using namespace std;
using namespace sdsl;

int main()
{
    csa_bitcompressed<> csa;
    construct_im(csa, "abracadabra", 1);
    cout << "csa.size(): " << csa.size() << endl;
    cout << "csa.sigma : " << csa.sigma << endl;
    cout << csa << endl;
    cout << extract(csa, 0, csa.size()-1) << endl;
}
