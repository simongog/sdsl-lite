//! Program for debugging. Load vector v of type int_vector<> from
// file argv[1]. Then the program outputs for each integer i,
// which is read from stdin, the value v[i].
#include <sdsl/int_vector.hpp>
#include <iostream>
#include <string>

using namespace std;
using namespace sdsl;

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " file [lb] [rb]" << endl;
        cout << " Reads a serialized int_vector<> from disk." <<endl;
        cout << " Outputs elements in the range [lb..rb], if specified." << endl;
        cout << " Otherwise, reads indexes from stdin and outputs values. "<<endl;
    }
    int_vector<> v;
    load_from_file(v, argv[1]);
    cout<<"loaded int_vector<> containing "<<v.size()<<" "<<
        (int)v.width()<<"-bit integers"<<endl;
    if (argc>3) {
        size_t a=stoull(argv[2]);;
        size_t b=stoull(argv[3]);
        if (b >= v.size())
            b = v.size()-1;
        if (a > b)
            a = b;
        for (size_t i=a; i<=b; ++i) {
            cout<<"v["<<i<<"]="<<v[i]<<endl;
        }
    } else {
        cout << "Interactive mode." << endl;
        size_t i;
        while (cin>>i) {
            cout<<"v["<<i<<"]="<<v[i]<<endl;
        }
    }
}
