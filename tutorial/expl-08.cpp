#include <iostream>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace sdsl;

int main()
{
    bit_vector b = bit_vector(8000, 0);
    for (size_t i=0; i < b.size(); i+=100)
        b[i] = 1;
    sd_vector<> sdb(b);
    sd_vector<>::rank_1_type sdb_rank(&sdb);
    for (size_t i=0; i<=sdb.size(); i+= sdb.size()/4)
        cout << "(" << i << ", " << sdb_rank(i) << ") ";
    cout << endl;
    rrr_vector<> rrrb(b);
    rrr_vector<>::rank_1_type rrrb_rank(&rrrb);
    for (size_t i=0; i<=rrrb.size(); i+= rrrb.size()/4)
        cout << "(" << i << ", " << rrrb_rank(i) << ") ";
    cout << endl;
}
