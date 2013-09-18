#include <sdsl/suffix_arrays.hpp>
#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>

using namespace sdsl;
using namespace std;

int main(int argc, char** argv)
{
    if (argc <  2) {
        cout << "Usage " << argv[0] << " text_file [max_locations] [post_context] [pre_context]" << endl;
        cout << "    This program constructs a very compact FM-index" << endl;
        cout << "    which supports count, locate, and extract queries." << endl;
        cout << "    text_file      Original text file." << endl;
        cout << "    max_locations  Maximal number of location to report." <<endl;
        cout << "    post_context   Maximal length of the reported post-context." << endl;
        cout << "    pre_context    Maximal length of the pre-context." << endl;
        return 1;
    }
    size_t max_locations = 5;
    size_t post_context = 10;
    size_t pre_context = 10;
    if (argc >= 3) {
        max_locations = atoi(argv[2]);
    }
    if (argc >= 4) {
        post_context = atoi(argv[3]);
    }
    if (argc >= 5) {
        pre_context = atoi(argv[4]);
    }
    string index_suffix = ".fm9";
    string index_file   = string(argv[1])+index_suffix;
    csa_wt<wt_huff<rrr_vector<127> >, 512, 1024> fm_index;

    if (!load_from_file(fm_index, index_file)) {
        ifstream in(argv[1]);
        if (!in) {
            cout << "ERROR: File " << argv[1] << " does not exist. Exit." << endl;
            return 1;
        }
        cout << "No index "<<index_file<< " located. Building index now." << endl;
        construct(fm_index, argv[1], 1); // generate index
        store_to_file(fm_index, index_file); // save it
    }
    cout << "Index construction complete, index requires " << size_in_mega_bytes(fm_index) << " MiB." << endl;
    cout << "Input search terms and press Ctrl-D to exit." << endl;
    string prompt = "\e[0;32m>\e[0m ";
    cout << prompt;
    string query;
    while (getline(cin, query)) {
        size_t m  = query.size();
        size_t occs = sdsl::count(fm_index, query.begin(), query.end());
        cout << "# of occurrences: " << occs << endl;
        if (occs > 0) {
            cout << "Location and context of first occurrences: " << endl;
            auto locations = locate(fm_index, query.begin(), query.begin()+m);
            sort(locations.begin(), locations.end());
            for (size_t i = 0, pre_extract = pre_context, post_extract = post_context; i < min(occs, max_locations); ++i) {
                cout << setw(8) << locations[i] << ": ";
                if (pre_extract > locations[i]) {
                    pre_extract = locations[i];
                }
                if (locations[i]+m+ post_extract > fm_index.size()) {
                    post_extract = fm_index.size()-locations[i]-m;
                }
                auto s   = extract(fm_index, locations[i]-pre_extract, locations[i]+m+ post_extract-1);
                string pre = s.substr(0, pre_extract);
                s = s.substr(pre_extract);
                if (pre.find_last_of('\n') != string::npos) {
                    pre = pre.substr(pre.find_last_of('\n')+1);
                }
                cout << pre;
                cout << "\e[1;31m";
                cout << s.substr(0, m);
                cout << "\e[0m";
                string context = s.substr(m);
                cout << context.substr(0, context.find_first_of('\n')) << endl;
            }
        }
        cout << prompt;
    }
    cout << endl;
}

