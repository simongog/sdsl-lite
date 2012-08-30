#include <sdsl/csa_wt.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/rmq_support.hpp>
#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <queue>

using namespace sdsl;
using namespace std;

struct interval {
    size_t lb, rb; // left bound, right bound (both inclusive)
    size_t min_val, min_idx; // minimal value, index of the minimum

    interval(size_t lb, size_t rb, size_t min_val, size_t min_idx):
        lb(lb), rb(rb), min_val(min_val), min_idx(min_idx) {};

    interval(const interval& i):
        lb(i.lb), rb(i.rb), min_val(i.min_val), min_idx(i.min_idx) {};

    bool operator>(const interval& i)const {
        if (min_val != i.min_val) {
            return min_val > i.min_val;
        }
        if (min_idx != i.min_idx) {
            return min_idx > i.min_idx;
        }
        if (lb != i.lb) {
            return lb > i.lb;
        }
        return rb > i.rb;
    }
};

int main(int argc, char** argv)
{
    if (argc <  2) {
        cout << "Usage " << argv[0] << " text_file [max_locations] [post_context] [pre_context]" << endl;
        cout << "      This program constructs a very compact FM-index" << endl;
        cout << "      which supports count, locate, and extract queries." << endl;
        cout << "      text_file      Original text file." << endl;
        cout << "      max_locations  Maximal number of location to report." <<endl;
        cout << "      post_context   Maximal length of the reported post-context." << endl;
        cout << "      pre_context    Maximal length of the pre-context." << endl;
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
    typedef csa_wt<wt_huff<rrr_vector<255> >, 512, 1024> fm_index_type;
    typedef fm_index_type::size_type size_type;
    fm_index_type fm_index;
    if (!util::load_from_file(fm_index, index_file.c_str())) {
		ifstream in(argv[1]);
		if( !in ){ cout << "ERROR: File " << argv[1] << " does not exist. Exit." << endl; return 1; }
		cout << "No index "<<index_file<< " located. Building index now." << endl;
        construct_csa(argv[1], fm_index); // generate index
        util::store_to_file(fm_index, index_file.c_str()); // save it
    }
    string rmq_suffix = ".rmq";
    string rmq_file   = string(argv[1])+rmq_suffix;
    rmq_succinct_sct<> rmq;
    csa_wt<wt_huff<>, 1,1<<20 > tmp_csa;
    if (!util::load_from_file(rmq, rmq_file.c_str())) {
		ifstream in(argv[1]);
		if( !in ){ cout << "ERROR: File " << argv[1] << " does not exist. Exit." << endl; return 1; }
		cout << "No index "<<rmq_file<< " located. Building index now." << endl;
        construct_csa(argv[1], tmp_csa);
        util::assign(rmq, rmq_succinct_sct<>(&tmp_csa));
        util::store_to_file(rmq, rmq_file.c_str());
    }
	size_t index_size = util::get_size_in_mega_bytes(fm_index) + util::get_size_in_mega_bytes(rmq);
	cout << "Index construction complete, index requires " << index_size << " MiB." << endl;
	cout << "Input search terms and press Ctrl-D to exit." << endl;
    string prompt = "\e[0;32m>\e[0m ";
    cout << prompt;
    string query;
    while (getline(cin, query)) {
        size_type m    = query.size();
        size_type lb = 0, rb=0; // left and right bound of the matching interval in the suffix array
        size_t occs = algorithm::backward_search(fm_index, 0, fm_index.size()-1, (const unsigned char*)query.c_str(), m, lb, rb);
        cout << "# of occcurences: " << occs << endl;
        if (occs > 0) {
            cout << "Location and context of first occurrences: " << endl;
            size_t min_idx = rmq(lb, rb);
            size_t min_val = fm_index[min_idx];
            priority_queue<interval, vector<interval>, greater<interval> > pq; // min-heap
            pq.push(interval(lb, rb, min_val, min_idx));
            for (size_t i = 0, pre_extract = pre_context, post_extract = post_context; i < min(occs, max_locations); ++i) {
                interval r = pq.top(); pq.pop();
                size_t location = r.min_val;
                if (r.min_idx > r.lb) {
                    min_idx = rmq(r.lb, r.min_idx-1);
                    pq.push(interval(r.lb, r.min_idx-1, fm_index[min_idx], min_idx));
                }
                if (r.min_idx < r.rb) {
                    min_idx = rmq(r.min_idx+1, r.rb);
                    pq.push(interval(r.min_idx+1, r.rb, fm_index[min_idx], min_idx));
                }
                cout << setw(8) << location << ": ";
                if (pre_extract > location) {
                    pre_extract = location;
                }
                if (location+m+ post_extract > fm_index.size()) {
                    post_extract = fm_index.size()-location-m;
                }
                string s   = algorithm::extract(fm_index, location-pre_extract, location+m+ post_extract-1);
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

