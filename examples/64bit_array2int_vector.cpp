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
    size_t x = util::file_size(argv[1]);
    const int BPI=8;
    cout<<"file size in bytes = "<<x<<endl;
    if (x % BPI != 0) {
        cout << "Error: x%"<<BPI<<" != 0"<<endl;
        return 1;
    }
//	(1) scan file and determine the largest value
    FILE* f = fopen(argv[1], "rb"); // open input file
    if (f == NULL) {
        cout<<"ERROR: could not open file "<<argv[1]<<endl;
        return 1;
    }
    const size_t BUF_SIZE = 6400000; // BUF_SIZE has to be a multiple of 64
    uint64_t* buf = (uint64_t*)malloc(BUF_SIZE*BPI);
    uint64_t max=0;
    for (size_t i=0, len=BUF_SIZE*BPI; i<x; i+=len) {
        len = BUF_SIZE*BPI;
        if (i+len > x) {
            len = x-i;
        }
        fread((char*)buf, 1, len, f);
        for (size_t j=0; j<len/BPI; ++j) {
            if (buf[j] > max)
                max = buf[j];
//			cout<<" "<<buf[j]<<endl;
        }
    }
    cout<<"Max value: "<<max<<endl;
    uint8_t width = bits::hi(max)+1;
    cout<<"width="<<(int)width<<endl;

//  (2) scan file, bit-compress values and write to outfile
    rewind(f); // reset file pointer
    string ofile = string(argv[1])+".int_vector";
    FILE* of = fopen(ofile.c_str(),"wb"); // open output file
    if (of == NULL) {
        cout<<"ERROR: could not open output file "<<argv[1]<<endl;
        return 1;
    }
    uint64_t bitlen=(x/BPI)*width;
    fwrite((char*)&bitlen, 8, 1, of);
    fwrite((char*)&width, 1, 1, of);
    int_vector<> v(BUF_SIZE, 0, width);
    for (size_t i=0, len=BUF_SIZE*BPI; i<x; i+=len) {
        len = BUF_SIZE*BPI;
        if (i+len > x) {
            len = x-i;
        }
        fread((char*)buf, 1, len, f);
        for (size_t j=0; j<len/BPI; ++j) {
            v[j] = buf[j];
        }
        fwrite((char*)v.data(), 8, ((len/BPI*width+63)/64), of);
    }
    free(buf);
    fclose(f);
    fclose(of);
    util::load_from_file(v, ofile.c_str());
    cout<<"v.size()="<<v.size()<<endl;
    cout<<"v[0]="<<v[0]<<endl;
    const bool do_check = false;
    if (do_check) {
        int_vector<> check;
        util::load_vector_from_file(check, argv[1], BPI);
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
