#include <sdsl/wavelet_trees.hpp>
#include <iostream>

using namespace std;
using namespace sdsl;

int main()
{
    wt_blcd<> wt;
    construct(wt, "expl-13.cpp", 1);
    for (size_t i=0; i < wt.size() and wt[i]!='\n'; ++i)
        cout << wt[i];
    cout << endl;
    cout << "number of lines  : " << wt.rank(wt.size(), '\n') << endl;
    cout << "first '=' in line: " << wt.rank(wt.select(1, '='),'\n')+1 << endl;
}
