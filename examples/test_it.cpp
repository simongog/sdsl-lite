#include <sdsl/wavelet_trees.hpp>
#include <iostream>
#include <bitset>

using namespace std;
using namespace sdsl;

int main()
{
    std::cout << bits::hi(2) << endl;
    std::cout << bits::hi(4) << endl;
    std::cout << bits::hi(8) << endl;
    std::cout << bits::hi(16) << endl;
    /*
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

/*
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