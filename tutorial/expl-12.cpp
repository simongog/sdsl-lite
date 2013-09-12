#include <iostream>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/vectors.hpp>

using namespace std;
using namespace sdsl;

int main()
{
    sd_vector<> sd_b = bit_vector {1,0,1,1,1,0,1,1,0,0,1,0,0,1};
    size_t ones = sd_vector<>::rank_1_type(&sd_b)(sd_b.size());
    sd_vector<>::select_1_type sdb_sel(&sd_b);

    cout << sd_b << endl;

    for (size_t i=1; i <= ones; ++i)
        cout << sdb_sel(i) << " ";
    cout << endl;
}


