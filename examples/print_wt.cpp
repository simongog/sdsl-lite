#include <string>
#include <vector>
#include <queue>
#include <array>
#include <sdsl/wavelet_trees.hpp>

using namespace sdsl;
using namespace std;

template<typename t_wt>
void visualize_wt_rec(const t_wt& wt, typename t_wt::node_type v, size_t level, vector<string>& out)
{
    if (!wt.is_leaf(v)) {
        if (out.size() < level+4) {
            out.push_back("");
            out.push_back("");
        }
        while (out[level+2].size() < out[level].size()) {
            out[level+2] += " ";
            out[level+3] += " ";
        }

        auto vs = wt.expand(v);
        if (!wt.empty(vs[0])) {
            visualize_wt_rec(wt, vs[0], level+2, out);
        }
        if (!wt.empty(vs[0]) and !wt.empty(vs[1])) {
            out[level+2] += " ";
            out[level+3] += " ";
        }
        if (!wt.empty(vs[1])) {
            visualize_wt_rec(wt, vs[1], level+2, out);
        }

        size_t begin = out[level].size();
        size_t end   = out[level+2].size();
        size_t size  = wt.size(v);
        size_t delta = (end-begin)-size;

        for (size_t i=0; i < delta/2; ++i) {
            out[level] += " ";
            out[level+1] += " ";
        }
        auto seq_vec = wt.seq(v);
        auto bit_it = wt.bit_vec(v).begin();
        for (auto it = seq_vec.begin(); it!=seq_vec.end(); ++it, ++bit_it) {
            out[level]   += *it;
            out[level+1] += *bit_it ? "1" : "0";
        }

        for (size_t i=0; i < (delta+1)/2; ++i) {
            out[level] += " ";
            out[level+1] += " ";
        }
    } else {
        auto seq = wt.seq(v);
        for (auto it = seq.begin(); it!=seq.end(); ++it) {
            out[level] += *it;
            out[level+1] += " ";
        }
    }
}

template<typename t_wt>
void visualize_wt(string s, string label)
{
    t_wt wt;
    construct_im(wt, s, 1);
    vector<string> vs(2,"");
    visualize_wt_rec(wt, wt.root(), 0, vs);
    cout << label << endl << endl;
    for (size_t i=0; i<vs.size(); ++i)
        cout<<vs[i]<<endl;
}

int main(int argc, char* argv[])
{
//    string s = "barbarabierbarbarbar";
    string s = "rhabarberbarbarabarbarbar";
    if (argc > 1) {
        s = argv[1];
    }
    cout << "T=" << s << endl;
    cout <<"--"<<endl;
    visualize_wt<wt_blcd<>>(s,"Balanced shape");
    visualize_wt<wt_huff<>>(s,"Huffman shape");
    visualize_wt<wt_hutu<>>(s,"Hu-Tucker shape");
    visualize_wt<wt_int<>>(s,"Balanced shape; fixed codeword length");
}
