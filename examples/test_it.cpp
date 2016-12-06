#include <sdsl/wavelet_trees.hpp>
#include <iostream>
#include <bitset>
#include <sdsl/k2_tree_comp_helper.hpp>
#include <sdsl/k2_tree_hybrid.hpp>

using namespace std;
using namespace sdsl;
using namespace k2_tree_ns;

template <typename t_x>
t_x createBitmask(t_x start, t_x end){
    t_x result = 0;
    for (auto i = start; i <= end ; ++i) {
        result |= (1 << i);
    }
    return result;
}

uint64_t createBitmask(int64_t start, int64_t end){
    uint64_t result = 0;
    for (auto i = start; i < end ; ++i) {
        result |= (1ULL<< i);
    }
    return result;
}

uint64_t test(const std::pair<uint, uint> &lhs) {
    const uint bitsToInterleaveForK2 = bits::hi(2) * 0;
    const uint bitsToInterleaveForKLeaves = bits::hi(4) * 1;

    //bitsOfMaximalValue might be < 8*max(sizeof(t_x),sizeof(t_y))
    const int bits = 32; //FIXME: only 32 bit for now

    auto rK1 = bitsToInterleaveForK2+bitsToInterleaveForKLeaves;
    auto lK1 = 2*rK1;

    auto lK2_f = bits - bitsToInterleaveForK2 - bitsToInterleaveForKLeaves;
    auto rK2_f = bits - bitsToInterleaveForK2;
    auto lK2   = 2*bitsToInterleaveForKLeaves;

    cout << "lK2_f " << lK2_f << endl;
    cout << "rK2_f " << rK2_f << endl;
    cout << "lK2 " << lK2 << endl;

    //set to one between 2*(bitsToInterleaveForK2+bitsToInterleaveForKLeaves) and 2*bitsOfMaximalValue
    uint64_t k_leaves_bitmask = createBitmask(0, 2*(bitsToInterleaveForKLeaves));
    uint64_t k_leaves_pre_bitmask = createBitmask(0, bitsToInterleaveForKLeaves);
/*
    std::cout << "k1_bitmask " << std::bitset<40>(k1_bitmask) << std::endl;
    std::cout << "k2_bitmask " << std::bitset<40>(k2_bitmask) << std::endl;
    std::cout << "kl_bitmask " << std::bitset<40>(k_leaves_bitmask) << std::endl;

    std::cout << "k1_pre_bitmask " << std::bitset<40>(k1_pre_bitmask) << std::endl;
    std::cout << "k2_pre_bitmask " << std::bitset<40>(k2_pre_bitmask) << std::endl;
    std::cout << "kl_pre_bitmask " << std::bitset<40>(k_leaves_pre_bitmask) << std::endl;
*/

    auto first = (interleave<2>::bits(lhs.first >> rK1, lhs.second >> rK1) << lK1);
    auto second = (interleave<2>::bits((lhs.first << lK2_f) >> rK2_f, (lhs.second << lK2_f) >> rK2_f) << lK2);
    auto third = (interleave<4>::bits(lhs.first & k_leaves_pre_bitmask, lhs.second & k_leaves_pre_bitmask) & k_leaves_bitmask);

    std::cout << "1: "<< bitset<8>(first) << endl;
    std::cout << "2: " << bitset<8>(second) << endl;
    std::cout << "3: " << bitset<8>(third) << endl;

    return first | second | third;




    /*uint64_t k1_bitmask = 17179852800;
    uint64_t k2_bitmask = 16320;
    uint k_leaves_bitmask = 63;
    uint k1_pre_bitmask = 130944;
    uint k2_pre_bitmask = 120;
    uint k_leaves_pre_bitmask = 63;

    auto first = interleave<4>::bits(lhs.first >> 7, lhs.second >> 7) << 14;
    auto second = interleave<2>::bits((lhs.first << 25) >> 29, (lhs.second << 25) >> 29) << 6;
    auto third = interleave<8>::bits(lhs.first & k_leaves_pre_bitmask, lhs.second & k_leaves_pre_bitmask);

    std::cout << "first "<< bitset<34>(first) << endl;
    std::cout << "second " << bitset<34>(second) << endl;
    std::cout << "third " << bitset<34>(third) << endl;

    first  = first & k1_bitmask;
    second = second & k2_bitmask;
    third  = third & k_leaves_bitmask;

    std::cout << "first "<< bitset<34>(first) << endl;
    std::cout << "second " << bitset<34>(second) << endl;
    std::cout << "third " << bitset<34>(third) << endl;

    return first | second | third;*/
}


//Parallel Bits Deposit
//x    HGFEDCBA
//mask 01100100
//res  0CB00A00
//x86_64 BMI2: PDEP
template <typename Integral>
Integral deposit_bits(Integral x, Integral mask) {
    Integral res = 0;
    for(Integral bb = 1; mask != 0; bb += bb) {
        if(x & bb) {
            res |= mask & (-mask);
        }
        mask &= (mask - 1);
    }
    return res;
}


template<unsigned long N>
struct bin {
    enum { value = (N%10)+2*bin<N/10>::value };
} ;

template<>
struct bin<0> {
    enum { value = 0 };
} ;


int main(int argc, char *argv[]) {

    std::string file_name = argv[1];
    typedef k2_tree_hybrid<4, 5, 2, 8, bit_vector, bit_vector> tested_type;
    tested_type k2tree;
    tested_type k2tree_correct;

    k2tree_correct.load_from_ladrabin(file_name);

    k2tree.load_from_ladrabin_construct_external(file_name, k2tree_correct);

    if (!(k2tree == k2tree_correct)){
        std::cout << "Trees differ" << std::endl;
    }

    /*uint mask = bin<1111000011001100101>::value;//1010101010...
    uint x = 0x000000FF;//1010101010...

    auto asd = deposit_bits(x, mask);
    std::cout << "Afet pdep: " << bitset<32>(asd) << std::endl;*/
    /*
    bit_vector tmp(4, 0);
    //bit_vector tmp1(32, 1);

    std::cout << tmp.size() << std::endl;
    std::copy(tmp1.begin(), tmp1.end(), tmp.begin());


    for (int i = 0; i < tmp.size(); i++){
        std::cout << tmp[i] << std::endl;
    }
*/
    /*
    std::cout << bitset<32>(64) << " " << bitset<32>(43) << std::endl;
    std::cout << bitset<32>(95) << " " << bitset<32>(155) << std::endl;

    std::cout << test(make_pair((uint)1,(uint)2)) << std::endl;
    std::cout << test(make_pair((uint)1,(uint)4)) << std::endl;



    int asd = createBitmask(5,20);
    std::cout << bitset<32>(asd) << endl;

    std::cout << bits::hi(2) << endl;
    std::cout << bits::hi(4) << endl;
    std::cout << bits::hi(8) << endl;
    std::cout << bits::hi(16) << endl;

    int_vector<> asd = {1,3,4,11,111};
    for (size_t i = 0; i < asd.size(); ++i) {
        std::cout << "in asd: " << asd.get_int(i*64) << std::endl;
    }

    std::string tmp_file = ram_file_name(util::to_string(util::pid()) + "_" + util::to_string(util::id()));
    int_vector_buffer<> dictionary_buffer(tmp_file, std::ios::out);
    for (size_t i = 0; i < asd.size(); ++i) {
        dictionary_buffer.push_back(asd[i]);
    }

    std::shared_ptr<int_vector<>> dictionary;
    dictionary_buffer.close();
    {
        int_vector<> tmp;
        load_from_file(tmp, tmp_file);
        std::cout << "Tmp size: " << tmp.size() << std::endl;
        for (size_t i = 0; i < tmp.size(); ++i) {
            std::cout << "in tmp: " << tmp[i] << std::endl;
        }
        dictionary.reset(new int_vector<>(tmp));
    }
    remove(tmp_file);

    std::cout << "Width" << std::to_string(dictionary->width()) << std::endl;
    std::cout << "Dicitonary size: " << dictionary->size() << std::endl;
    for (size_t i = 0; i < dictionary->size(); ++i) {
        std::cout << "in dict: " << dictionary->operator[](i) << std::endl;
    }*/

/*    wt_huff_int<> wt;
    int_vector<> vec = int_vector<>(9, 0, 64);
    util::set_to_id(vec);  // 0 1 2 3 4 5 6 7 8
    //{3,12,4,4,5,1,6,4,2};
    construct_im(wt, vec, 8);
    cout << wt << endl;


    wt_huff_int<> wt;
    int_vector<> vec = int_vector<>(9, 0);
    util::set_to_id(vec);  // 0 1 2 3 4 5 6 7 8
    util::bit_compress(vec);
    //{3,12,4,4,5,1,6,4,2};
    construct_im(wt, vec);
    cout << wt << endl;


    for (size_t i=0; i < wt.size() and wt[i]!='\n'; ++i)
        cout << std::to_string(wt[i]) << "\t";
    cout << endl;
*/

/*
    int_vector<64> x_vec(6, 0);
    util::set_to_id(x_vec);
    cout << x_vec << endl;  // 0 1 2 3 4 5
    wt_huff_int<> wt;
    construct_im(wt, x_vec.raw, 8);
    cout << "wt.sigma : " << wt.sigma << endl;
    cout << "wt.size : " << wt.size() << endl;
    cout << wt << endl; // 0 1 2 3 4 5*/

}