#include <sdsl/int_vector.hpp>
#include <sdsl/suffix_trees.hpp>
#include <iostream>

using namespace std;
using namespace sdsl;

int main(int argc, char* argv[])
{

    bit_vector b(10000000, 0);
    b[8] = 1;
    rank_support_v<> rb(&b);

    cout<<rb(8)<<endl;
    cout<<rb(9)<<endl;

    cout<< "size of b in MB: " << size_in_mega_bytes(b)<< endl;
    cout<< "size of rb in MB: " << size_in_mega_bytes(rb)<< endl;

    rrr_vector<127> rrrb(b);
    rrr_vector<127>::rank_1_type rank_rrrb(&rrrb);
    cout<<rank_rrrb(8)<<endl;
    cout<<rank_rrrb(9)<<endl;

    cout<< "size of rrrb in MB: " << size_in_mega_bytes(rrrb)<< endl;
    cout<< "size of rank_rrrb in MB: " << size_in_mega_bytes(rank_rrrb)<< endl;


    rrr_vector<127>::select_1_type select_rrrb(&rrrb);
    cout<<"position of first one in b: "<<select_rrrb(1)<<endl;

    bit_vector x;
    util::assign(x, bit_vector(10000000,1));

    int_vector<> v(100, 5, 7);

    cout<<"v[5]="<<v[5]<<endl;
    v[5]=120;
    cout<<"v[5]="<<v[5]<<endl;


    int_vector<32> w(100, 4);

    write_structure<JSON_FORMAT>(rrrb, cout);
    cout<<endl;

    typedef cst_sada<> tCst;
    tCst cst;

    construct(cst, argv[1], 1);

    for (tCst::const_iterator it = cst.begin(); it!=cst.end(); ++it) {
        if (it.visit() == 1) {
            cout << cst.depth(*it) << "-["<< cst.lb(*it)<<","<<cst.rb(*it)<<"]" << endl;

        }
    }
}
