/*
 * Transforms a sequence of 64-bit integers into a bit-compressed integer vector.
 * The first command line parameter argv[1] specifies the file, which contains the sequence
 * of integers.
 * The bit-compressed integer vector is stored in a file called `argv[1].int_vector`.
 */
#include <sdsl/util.hpp>
#include <sdsl/int_vector.hpp>
#include <iostream>
#include <string>

using namespace sdsl;
using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " int_file" << endl;
        return 1;
    }

    size_t x = util::file_size(argv[1]);
    const int BPI=8;
    cout<<"file size in bytes = "<<x<<endl;
    if (x % BPI != 0) {
        cout << "Error: x%"<<BPI<<" != 0"<<endl;
        return 1;
    }
//	(1) scan file and determine the largest value
    int_vector_buffer<64> ivb(string(argv[1]), std::ios::in, 6400000, 64, true);
    if (!ivb.is_open()) {
        cout<<"ERROR: could not open file "<<argv[1]<<endl;
        return 1;
    }
    uint64_t max=0;
    for (size_t i=0; i<ivb.size(); ++i) {
        if (ivb[i] > max)
            max = ivb[i];
    }
    cout<<"Max value: "<<max<<endl;
    uint8_t width = bits::hi(max)+1;
    cout<<"width="<<(int)width<<endl;

//  (2) scan file, bit-compress values and write to outfile
    string ofile = string(argv[1])+".int_vector";
    int_vector_buffer<> ivb2(ofile, std::ios::out, 6400000, width);
    if (!ivb2.is_open()) {
        cout<<"ERROR: could not open output file "<<argv[1]<<endl;
        return 1;
    }
    for (size_t i=0; i<ivb.size(); ++i) {
        ivb2[i] = ivb[i];
    }
    ivb.close();
    ivb2.close();

    int_vector<> v;
    load_from_file(v, ofile);
    cout<<"v.size()="<<v.size()<<endl;
    cout<<"v[0]="<<v[0]<<endl;
    const bool do_check = false;
    if (do_check) {
        int_vector<> check;
        load_vector_from_file(check, argv[1], BPI);
        if (check.size() != v.size()) {
            cout<<"Vectors differ in size: "<<check.size()<<"!="<<v.size()<<endl;
            return 1;
        }
        for (size_t i=0; i<check.size(); ++i) {
            if (check[i] != v[i]) {
                cout<<"vectors differ in position "<<i<<": "<<check[i]<<"!="<<v[i]<<endl;
            }
        }
    }
}
