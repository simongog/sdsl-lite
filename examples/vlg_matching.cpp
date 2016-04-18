#include "../include/sdsl/vlg_index.hpp"
//#include <sdsl/vlg_index.hpp>

using namespace sdsl;
using namespace std;

template<typename index>
void dump_query_results(index idx, string qry)
{
    cout << endl;
    cout << "count(" << qry << ")=" << count(idx, qry) << endl;
    cout << "locate(" << qry << ")=" << endl;
    auto res = locate(idx, qry);
    size_t occ = 1;
    for (auto it = res.begin(); it != res.end(); ++it, ++occ) {
        auto pos = *it;
        cout << "  " << occ << ". occ starting at position " << pos << endl;
        cout << "     Subpattern positions:";
        for (size_t i=0; i<it.size(); ++i) {
            cout << " " << it[i];
        }
        cout << endl;
    }
}

int main(int argc, char* argv[])
{
    // index with byte alphabet (default)
    {
        vlg_index<> idx;
        
        if (argc <= 1) {
            construct_im(idx, "abracadabrasimsalabim", 1);
        } else {
            string filename(argv[1]);
            construct(idx, filename, 1);
            ofstream out(filename + ".info.html");
            write_structure<HTML_FORMAT>(idx, out);
        }
        
        dump_query_results(idx, "ac.{2,5}?a.{4,8}?b");
    }
    
    // index with int alphabet
    {
        vlg_index<int_alphabet_tag> idx;
        
        if (argc <= 1) {
            construct_im(idx, int_vector<>({97, 98, 114, 97, 99, 97, 100, 97, 98, 114, 97, 115, 105, 109, 115, 97, 108, 97, 98, 105, 109}));
        } else {
            string filename(argv[1]);
            construct(idx, filename, 0);
            ofstream out(filename + ".info.html");
            write_structure<HTML_FORMAT>(idx, out);
        }
        
        dump_query_results(idx, "97 99 .{2,5}? 97 .{4,8}? 98");
    }
}
