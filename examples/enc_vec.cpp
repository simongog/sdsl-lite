#include <iostream>
#include <sdsl/vectors.hpp>

using namespace std;
using namespace sdsl;

int main()
{
    constexpr uint64_t d = 3;
    vector<uint64_t> v = {1,2,3,4,6,8,10,11,14,16,17,18,19,20};
//    vector<uint64_t> v = {1,2,3,4,6,8,10,11,14,16,17,18,19,20};
    {
        enc_vector<coder::elias_delta, d> ev(v);
        for (size_t i=0; i<v.size(); ++i) {
            cout << v[i] << " = " << ev[i] << endl;
        }
    }

    int_vector<> vv(v.size());
    for (size_t i=0; i<v.size(); ++i)
        vv[i] = v[i];
    util::bit_compress(vv);
    store_to_file(vv, "v.sdsl");
    {
        int_vector_buffer<> ivb("v.sdsl");
        enc_vector<coder::elias_delta, d> ev(ivb);
        for (size_t i=0; i<v.size(); ++i) {
            cout << v[i] << " = " << ev[i] << endl;
        }
    }
    sdsl::remove("v.sdsl");
    /*
        cout << "--" << endl;
        uint64_t store[d] = {0};
        for (size_t i=0; i*d < ev.size(); ++i) {
            ev.get_inter_sampled_values(i, store);
            cout<<" "<<ev.sample(i)<<endl;
            for (uint64_t j=0; j<d and i*d+j+1 < ev.size(); ++j){
                cout << ev.sample(i)+store[j] << " = " << ev[i*d+j] << "   " << v[j]<<endl;
            }
            cout <<"-"<<endl;
        }
        */
}
