#include <sdsl/bit_vectors.hpp>
#include <random>
#include <iostream>
#include <chrono>

using namespace sdsl;
using namespace std;

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;


template<class t_vec>
uint64_t test_inv_random_access(const t_vec& v, const int_vector<64>& rands, uint64_t mask, uint64_t times=100000000)
{
    uint64_t cnt=0;
    for (uint64_t i=0; i<times; ++i) {
        cnt += v(rands[ i&mask ]);
    }
    return cnt;
}



//int main(int argc, char* argv[]){
int main()
{
    auto start = timer::now();
    bool default_value = 0; //ID[ID.length()-1]-'0';
    bit_vector bv = bit_vector(800000000, default_value);

    std::mt19937_64 rng;
    std::uniform_int_distribution<uint64_t> distribution(0, bv.size()-1);
    auto dice = bind(distribution, rng);
    // populate vectors with some other bits
    for (uint64_t i=0; i < bv.size()/25; ++i) {
        uint64_t x = dice();
        bv[x] = !default_value;
    }
    auto stop = timer::now();
    cout << "initialization in (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
    cout << "size in MiB: " << size_in_mega_bytes(bv) << endl;

    start = timer::now();
    sd_vector<> bv_sd(bv);
    stop = timer::now();
    cout << "sd_construction in (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
    {
        bit_vector().swap(bv);
    }
    cout << "size in MiB: " << size_in_mega_bytes(bv_sd) << endl;
    cout << "wl = " << (size_t) bv_sd.wl << endl;
    cout << "n = " << bv_sd.size() << endl;
    cout << "2*m = " << bv_sd.high.size()<<endl;
    cout <<"n/m=" << (2.0*bv_sd.size())/bv_sd.high.size()<<endl;

    auto zeros = sd_vector<>::rank_0_type(&bv_sd)(bv_sd.size());
    auto ones = bv_sd.size()-zeros;
    cout << "zeros = "<< zeros << endl;
    {
        uint64_t mask = 0;
        auto rands = util::rnd_positions<int_vector<64>>(20, mask, zeros, 17);
        for (uint64_t i=0; i<rands.size(); ++i) rands[i] = rands[i]+1;
        sd_vector<>::select_0_type select0(&bv_sd);
        const uint64_t reps = 10000000;
        start = timer::now();
        auto check = test_inv_random_access(select0, rands, mask, reps);
        stop = timer::now();

        cout << "# select0_time = " << duration_cast<nanoseconds>(stop-start).count()/(double)reps << endl;
        cout << "# select_check = " << check << endl;
        cout << "# size_in_mega_bytes(bv_sd) = " << size_in_mega_bytes(bv_sd) << endl;
        cout << "# size_in_mega_bytes(select0) = " << size_in_mega_bytes(select0) << endl;
    }
    {
        uint64_t mask = 0;
        auto rands = util::rnd_positions<int_vector<64>>(20, mask, zeros, 17);
        for (uint64_t i=0; i<rands.size(); ++i) rands[i] = rands[i]+1;
        select_0_support_sd<sd_vector<>> select0(&bv_sd);
        const uint64_t reps = 10000000;
        start = timer::now();
        auto check = test_inv_random_access(select0, rands, mask, reps);
        stop = timer::now();

        cout << "# select0_time = " << duration_cast<nanoseconds>(stop-start).count()/(double)reps << endl;
        cout << "# select_check = " << check << endl;
        cout << "# size_in_mega_bytes(bv_sd) = " << size_in_mega_bytes(bv_sd) << endl;
        cout << "# size_in_mega_bytes(select0) = " << size_in_mega_bytes(select0) << endl;
    }
    {
        uint64_t mask = 0;
        auto rands = util::rnd_positions<int_vector<64>>(20, mask, ones, 17);
        for (uint64_t i=0; i<rands.size(); ++i) rands[i] = rands[i]+1;
        sd_vector<>::select_1_type select1(&bv_sd);
        const uint64_t reps = 10000000;
        start = timer::now();
        auto check = test_inv_random_access(select1, rands, mask, reps);
        stop = timer::now();

        cout << "# select1_time = " << duration_cast<nanoseconds>(stop-start).count()/(double)reps << endl;
        cout << "# select_check = " << check << endl;
    }
    {
        uint64_t mask = 0;
        auto rands = util::rnd_positions<int_vector<64>>(20, mask, bv_sd.size(), 17);
        cout<<"done"<<endl;
        cout<<(uint64_t)&(bv_sd.high_1_select)<<endl;
        cout<<(uint64_t)&(bv_sd.high_0_select)<<endl;
        sd_vector<>::rank_1_type rank1(&bv_sd);
        cout<<"done"<<endl;
        const uint64_t reps = 10000000;
//        for(size_t i=0; i<bv_sd.size();++i){
//            cout << "i="<<i<<" rank1("<<i<<")="<<rank1(i)<<endl;
//        }
        start = timer::now();
        auto check = test_inv_random_access(rank1, rands, mask, reps);
        stop = timer::now();

        cout << "# rank1_time = " << duration_cast<nanoseconds>(stop-start).count()/(double)reps << endl;
        cout << "# select_check = " << check << endl;
    }
}
