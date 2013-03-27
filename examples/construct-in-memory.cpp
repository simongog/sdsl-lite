#include <sdsl/suffix_trees.hpp>
#include <string>

using namespace sdsl;
using namespace std;

int main()
{
    string file = "@testinput.iv8";
    string v = "abracadabradabarab";
    cout<<"v.size()="<<v.size()<<endl;
    cout<<v<<endl;
    store_to_file((const char*)v.c_str(), file);
    cout<<"---------"<<endl;
    cout << util::file_size(file) << endl;
    cout<<"--------- construct wt_huff<> ----------"<<endl;
    {
        wt_huff<> wt;
        construct(wt, file, 1);
        cout << "wt.size()="<<wt.size()<<endl;
        cout << util::file_size(file) << endl;
        for (size_t i=0; i < wt.size(); ++i) {
            cout << wt[i] << " ";
        }
        cout << endl;
    }
    cout<<"---------"<<endl;
    cout<<"--------- construct csa_wt<> ----------"<<endl;
    {
        csa_wt<> csa;
        construct(csa, file, 1);
        cout << "csa.size()="<<csa.size()<<endl;
        for (size_t i=0; i < csa.size(); ++i) {
            cout << csa[i] << " ";
        }
        cout << endl;
    }
    cout<<"---------"<<endl;

    cout<<"--------- construct csa_sct3<> ----------"<<endl;
    {
        typedef cst_sada<> cst_type;
        cst_type cst;
        construct(cst, file, 1);
        cout << "cst.size()="<<cst.size()<<endl;
        for (cst_type::const_iterator it=cst.begin(); it!=cst.end(); ++it) {
            cout << cst.depth(*it) << "-[" << cst.lb(*it) << "," << cst.rb(*it) << "]" << endl;
        }
        cout << endl;
    }
    cout<<"---------"<<endl;

}
