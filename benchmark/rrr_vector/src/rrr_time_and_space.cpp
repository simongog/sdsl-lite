#include <iostream>
#include <fstream>
#include <string>
#include <sdsl/rrr_vector.hpp>

using namespace std;
using namespace sdsl;

#ifndef BLOCK_SIZE
#define 31 BLOCK_SIZE
#endif

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

int main(int argc, char* argv[])
{
    if (argc < 3) {
        cout << "Usage: " << argv[0] << " bit_vector_file k" << endl;
        cout << " generates a rrr_vector<" << BLOCK_SIZE << "> with sample rate k" << endl;
        cout << " for the bitvector stored in bit_vector_file and run a benchmark" << endl;
        return 1;
    }

    typedef rrr_vector<BLOCK_SIZE> rrr_vec_type;
    typedef rrr_vec_type::select_1_type rrr_select_type;
    typedef rrr_vec_type::rank_1_type rrr_rank_type;
    bit_vector bv;
    if (load_from_file(bv, argv[1])) {
        cout << "# plain_size = " << size_in_bytes(bv) << endl;
        uint16_t k = atoi(argv[2]);
        auto start = timer::now();
        rrr_vec_type rrr_vector(bv, k);
        util::clear(bv);
        rrr_select_type rrr_sel(&rrr_vector);
        rrr_rank_type   rrr_rank(&rrr_vector);
        auto stop = timer::now();
        cout << "# construct_time = " << duration_cast<milliseconds>(stop-start).count() << endl;
        rrr_vec_type::size_type args = rrr_rank(rrr_vector.size());
        cout << "# rrr_vector.size() = " << rrr_vector.size() << endl;
        cout << "# args = " << args << endl;
        cout << "# file_name = "  << argv[1] << endl;
        cout << "# block_size = " << BLOCK_SIZE << endl;
        cout << "# sample_rate = "<< k << endl;
        cout << "# rrr_size = "   << size_in_bytes(rrr_vector) << endl;
        cout << "# bt_size = "    << size_in_bytes(rrr_vector.bt) << endl;
        cout << "# btnr_size = "  << size_in_bytes(rrr_vector.btnr) << endl;
        const uint64_t reps = 10000000;
        uint64_t mask = 0;
        int_vector<64> rands = util::rnd_positions<int_vector<64>>(20, mask, rrr_vector.size(), 17);
        start = timer::now();
        test_random_access(rrr_vector, rands, mask, reps);
        stop = timer::now();
        cout << "# access_time = " << duration_cast<nanoseconds>(stop-start).count()/(double)reps << endl;
        rands = util::rnd_positions<int_vector<64>>(20, mask, rrr_vector.size()+1, 17);
        start = timer::now();
        test_inv_random_access(rrr_rank, rands, mask, reps);
        stop = timer::now();
        cout << "# rank_time = " << duration_cast<nanoseconds>(stop-start).count()/(double)reps << endl;
        rands = util::rnd_positions<int_vector<64>>(20, mask, args, 17);
        for (uint64_t i=0; i<rands.size(); ++i) rands[i] = rands[i]+1;
        stop = timer::now();
        test_inv_random_access(rrr_sel, rands, mask, reps);
        stop = timer::now();
        cout << "# select_time = " << duration_cast<nanoseconds>(stop-start).count()/(double)reps << endl;
    }
}
