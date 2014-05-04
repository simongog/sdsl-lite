#include "sdsl/int_vector.hpp"
#include <string>
#include <cstdlib>
#include <random>

using namespace std;
using namespace sdsl;

int main(int argc, char* argv[])
{
    if (argc < 3) {
        cout << "Usage: " << argv[0] << " FILE BIT_VECTOR_ID" << endl;
        cout << " Generator program for bit_vectors." << endl;
        cout << " (1) The corresponding bit_vector for BIT_VECTOR_ID is generated." << endl;
        cout << " (2) The bit_vector is stored to FILE." << endl;
        return 1;
    }
    string ID = argv[2];
    bit_vector bv;
    if ("CRAFTED-32" == ID) {
        bv = bit_vector(32, 0);
        bv[1] = bv[4] = bv[7] = bv[18] = bv[24] = bv[26] = bv[30] = bv[31] = 1;
    } else if ("CRAFTED-SPARSE-0" == ID or "CRAFTED-SPARSE-1" == ID) {
        bool default_value = ID[ID.length()-1]-'0';
        bv = bit_vector(1000000, default_value);
        std::mt19937_64 rng;
        std::uniform_int_distribution<uint64_t> distribution(0, bv.size()-1);
        auto dice = bind(distribution, rng);
        // populate vectors with some other bits
        for (uint64_t i=0; i < bv.size()/1000; ++i) {
            uint64_t x = dice();
            bv[x] = !default_value;
        }
    } else if ("CRAFTED-BLOCK-0" == ID or "CRAFTED-BLOCK-1" == ID) {
        bool default_value = ID[ID.length()-1]-'0';
        std::mt19937_64 rng;
        std::uniform_int_distribution<uint64_t> distribution1(0, bv.size()-1);
        std::uniform_int_distribution<uint64_t> distribution2(0, 999);
        auto dice1 = bind(distribution1, rng);
        auto dice2 = bind(distribution2, rng);
        bv = bit_vector(1000000, default_value);
        // populate vectors with some blocks of other bits
        for (uint64_t i=0; i< bv.size()/1000; ++i) {
            uint64_t x = dice1();
            uint64_t len = dice2();
            for (uint64_t j=x; j<x+len and j<bv.size(); ++j) {
                bv[j] = 1-default_value;
            }
        }
    } else if ("CRAFTED-MAT-SELECT") {
        // Matthias Petri's test
        srand(4711);
        uint64_t ones = 4030+(rand()%80);
        bv = bit_vector(1000000);
        for (uint64_t i=0; i<ones; i++) {
            uint64_t rnd = rand();
            while (bv[rnd%bv.size()]==1) {
                rnd = rand();
            }
            bv[rnd%bv.size()] = 1;
        }
    } else {
        cerr << "ID not found. No bitvector generated." << endl;
        return 1;
    }
    store_to_file(bv, argv[1]);
}
