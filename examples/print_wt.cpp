#include <string>
#include <queue>
#include <array>
#include <sdsl/wavelet_trees.hpp>

using namespace sdsl;
using namespace std;

template<typename t_wt>
void console_print_wt(const t_wt& wt)
{
    auto v = wt.root();
    array<queue<decltype(v)>,2> q;
    uint64_t level = 0;
    q[level%2].push(v);
    while (!q[level%2].empty()) {
        while (!q[level%2].empty()) {
            v = q[level%2].front();
            q[level%2].pop();
            if (!wt.empty(v)) {
                if (!wt.is_leaf(v)) {
                    for (auto it = wt.begin(v); it!=wt.end(v); ++it) {
                        cout<<*it;
                    }
                    cout<<" ";
                    auto vs = wt.expand(v);
                    q[(level+1)%2].push(get<0>(vs));
                    q[(level+1)%2].push(get<1>(vs));
                } else {
                    for (size_t i=0; i<wt.size(v); ++i) {
                        cout<<wt.sym(v);
                    }
                    cout<<" ";
                }
            }
        }
        ++level;
        cout << endl;
    }
}

int main(int argc, char* argv[])
{
    string s = "barbarabierbarbarbar";
    if (argc > 1) {
        s = argv[1];
    }
    cout << "T=" << s << endl;
    wm_int<> wt;
//    wt_blcd<bit_vector, bit_vector::rank_1_type, bit_vector::select_1_type, bit_vector::select_0_type, int_tree<>> wt;
//    wt_hutu<> wt;
//    wt_huff<> wt;
//    wt_blcd<> wt;
//    wt_int<> wt;
    construct_im(wt, s, 1);
    cout <<"--"<<endl;
    cout << "T=" << wt << endl;

    console_print_wt(wt);
}
