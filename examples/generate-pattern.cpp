#include <sdsl/suffix_trees.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

typedef cst_sct3<csa_wt<wt_int<>>> cst_type;

int main(int argc, char** argv)
{
    if (argc <= 3) {
        cout << "Usage: " << argv[0] << " file count length [charset]" << endl;
        cout << " - file:       The data to extract the patterns from" << endl;
        cout << " - count:      Number of patterns to extract" << endl;
        cout << " - length:     Minimum length of patterns to extract" << endl;
        cout << " - human readable text: 0 or 1 (default)" << endl;
        cout << " - charset (only applies if human readable):" << endl;
        cout << "   - 0 = alpha-numerical only (default)" << endl;
        cout << "   - 1 = alpha-numerical with whitespace" << endl;
        cout << endl;
        cout << " Uses a CST to extract most common patterns/phrases matching provided criteria." << endl;
        cout << endl;
        return 1;
    }

    // parse args
    string file(argv[1]);
    size_t count = atoi(argv[2]);
    size_t length = atoi(argv[3]);
    bool human_readable = argc <= 4 ? true : (atoi(argv[4]) == 1);
    int charset = argc <= 5 ? 0 : atoi(argv[5]);

    // CHARSET
    auto is_special_char = [](char c) { return isspace(c) || ispunct(c); };

    std::function<bool(char c)> charset_predicates[] = {
        [](char c) { return isalnum(c) || c == '_'; },
        [](char c) { return isalnum(c) || c == '_' || c == ' '; }
    };
    auto charset_predicate = charset_predicates[charset];
    auto charset_string_predicate = [&](cst_type::string_type s) { return !human_readable || all_of(s.begin(), s.end(), charset_predicate); };

    // load (cached) CST
    string cst_file = file + ".sdsl";
    cst_type cst;
    if (!load_from_file(cst, cst_file)) {
        construct(cst, file, 0);
        store_to_file(cst, cst_file);
    }

    // traverse
    using entry_type = pair<size_t, cst_type::string_type>;
    priority_queue<entry_type, vector<entry_type>, greater<entry_type>> found;
    auto min_occ = [&]() { return found.empty() ? 0 : found.top().first; };
    for (auto it=cst.begin(); it!=cst.end(); ++it) {
        if (it.visit() == 1) {
            auto v = *it;
            auto depth = cst.depth(v);
            if (depth == 0)
                continue;

            auto curr_min_occ = min_occ();
            auto curr_occ = cst.size(v);
            if (curr_occ <= curr_min_occ && found.size() == count) {
                it.skip_subtree();
                continue;
            }

            auto begin = cst.csa[cst.lb(v)];
            cst_type::string_type phrase = extract(cst.csa, begin, begin + min((size_t)depth, (size_t)length) - 1);

            // check charset
            if (!charset_string_predicate(phrase)) {
                it.skip_subtree();
                continue;
            }

            // check length
            if (depth >= length) {
                // check for word boundary (approximation)
                /* if (human_readable
                     && !is_special_char(cst.csa.L[cst.rb(v)])
                     && !is_special_char(cst.csa.L[cst.lb(v)])) {
                     it.skip_subtree();
                     continue;
                 }*/

                found.emplace(curr_occ, phrase);
                if (found.size() > count)
                    found.pop();

                // report status
                cerr << "Found " << found.size() << " phrases with min. #occurrences of " << found.top().first << " so far" << endl;

                it.skip_subtree();
            }
        }
    }

    // fill
    while (found.size() < count)
        found.emplace(0, cst_type::string_type(length, 'x'));

    // output
    while (!found.empty()) {
        auto symbols = found.top().second;
        if (human_readable)
            cout << string(symbols.begin(), symbols.end()) << endl;
        else
            cout << symbols << endl;
        found.pop();
    }
}
