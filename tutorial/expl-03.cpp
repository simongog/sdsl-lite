#include <iostream>
#include <sdsl/vectors.hpp>
#include <sdsl/coder_elias_gamma.hpp>
#include <sdsl/coder_comma.hpp>
#include <sdsl/coder_fibonacci.hpp>

using namespace std;
using namespace sdsl;

int main()
{
    int_vector<> v(10*(1<<20), 0);
    v[0] = 1ULL<<63;
    //util::bit_compress(v);
    cout << size_in_mega_bytes(v) << endl;
    vlc_vector<coder::fibonacci> vv(v);
    cout << size_in_mega_bytes(vv) << endl;
	cout << "Percentage: " << size_in_mega_bytes(vv) / size_in_mega_bytes(v) * 100 << endl;
}
