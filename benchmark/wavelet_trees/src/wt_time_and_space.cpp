#include<iostream>
#include<fstream>
#include<string>
#include<sdsl/wavelet_trees.hpp>
#include<sdsl/wt_helper.hpp>

using namespace std;
using namespace sdsl;
using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

typedef int_vector<>::size_type size_type;
typedef WT_TYPE::value_type value_type;


// test access
template<class t_wt>
uint64_t test_access(const t_wt& wt, const vector<size_type>& is, uint64_t mask, uint64_t times=100000000)
{
    uint64_t cnt=0;
    for (uint64_t i=0; i<times; ++i) {
        cnt += wt[is[i&mask]];
    }
    return cnt;
}

// test rank
template<class t_wt>
uint64_t test_rank(const t_wt& wt, const vector<size_type>& is, const vector<value_type>& cs, uint64_t mask, uint64_t times=100000000)
{
    uint64_t cnt=0;
    for (uint64_t i=0; i<times; ++i) {
        cnt += wt.rank(is[i&mask], cs[i&mask]);
    }
    return cnt;
}

// test inverse_select
template<class t_wt>
uint64_t test_inverse_select(const t_wt& wt, const vector<size_type>& is, uint64_t mask, uint64_t times=100000000)
{
    uint64_t cnt=0;
    for (uint64_t i=0; i<times; ++i) {
        cnt += wt.inverse_select(is[i&mask]).first;
    }
    return cnt;
}

// test interval_symbols
template<class t_wt>
uint64_t
test_interval_symbols(typename enable_if<!(has_node_type<t_wt>::value),
                      t_wt>::type&, const vector<size_type>&, const vector<size_type>&, size_type&, uint64_t, uint64_t)
{
    return 0; // interval_symbols not implemented
}

template<class t_wt>
uint64_t
test_interval_symbols(typename enable_if<has_node_type<t_wt>::value,
                      t_wt>::type& wt, const vector<size_type>& is, const vector<size_type>& js, size_type& k, uint64_t mask, uint64_t times=100000000)
{
    vector<value_type> tmp(wt.sigma);
    vector<size_type> tmp2(wt.sigma);
    uint64_t cnt=0;
    for (uint64_t i=0; i<times; ++i) {
        interval_symbols(wt, is[i&mask], js[i&mask], k, tmp, tmp2, tmp2);
        cnt += k;
    }
    return cnt;
}

// test lex_count
template<class t_wt>
uint64_t
test_lex_count(typename enable_if<!(t_wt::lex_ordered), t_wt>::type&, const vector<size_type>&, const vector<size_type>&, const vector<value_type>&, uint64_t, uint64_t)
{
    return 0; // lex_count not implemented
}

template<class t_wt>
uint64_t
test_lex_count(typename enable_if<t_wt::lex_ordered, t_wt>::type& wt, const vector<size_type>& is, const vector<size_type>& js, const vector<value_type>& cs, uint64_t mask, uint64_t times=100000000)
{
    uint64_t cnt=0;
    for (uint64_t i=0; i<times; ++i) {
        cnt += get<0>(wt.lex_count(is[i&mask], js[i&mask], cs[i&mask]));
    }
    return cnt;
}

// test lex_smaller_count
template<class t_wt>
uint64_t
test_lex_smaller_count(typename enable_if<!(t_wt::lex_ordered), t_wt>::type&, const vector<size_type>&, const vector<value_type>&, uint64_t, uint64_t)
{
    return 0; // lex_smaller_count not implemented
}

template<class t_wt>
uint64_t
test_lex_smaller_count(typename enable_if<t_wt::lex_ordered, t_wt>::type& wt, const vector<size_type>& is, const vector<value_type>& cs, uint64_t mask, uint64_t times=100000000)
{
    uint64_t cnt=0;
    for (uint64_t i=0; i<times; ++i) {
        cnt += get<0>(wt.lex_smaller_count(is[i&mask], cs[i&mask]));
    }
    return cnt;
}

// test select
template<class t_wt>
uint64_t test_select(const t_wt& wt, const vector<size_type>& is, const vector<value_type>& cs, uint64_t mask, uint64_t times=100000000)
{
    uint64_t cnt=0;
    for (uint64_t i=0; i<times; ++i) {
        cnt += wt.select(is[i&mask], cs[i&mask]);
    }
    return cnt;
}

// generate benchmark input
template<class t_iv>
void random_cs(const t_iv& iv, vector<value_type>& cs)
{
    std::mt19937_64 rng;
    std::uniform_int_distribution<uint64_t> distribution(0, iv.size()-1);
    auto dice = bind(distribution, rng);
    for (uint64_t l = 0; l<cs.size(); ++l) {
        cs[l]=iv[dice()];
    }
}

template<class t_iv>
void random_is_js(const t_iv& iv, vector<size_type>& is, vector<size_type>& js)
{

    std::mt19937_64 rng;
    std::uniform_int_distribution<uint64_t> distribution(0, iv.size()-1);
    auto dice = bind(distribution, rng);
    for (uint64_t l = 0; l<is.size(); ++l) {
        is[l]=dice();
        js[l]=min(is[l]+dice(), iv.size());
    }
}

template<class t_iv>
void prepare_for_select(const t_iv& iv, vector<value_type>& cs, vector<size_type>& is)
{

    uint64_t sigma = 0;
    for (uint64_t i=0; i<iv.size(); ++i) {
        if (sigma < iv[i]) sigma = iv[i];
    }
    ++sigma;

    vector<size_type> symbols(sigma, 0);

    for (uint64_t i=0; i<iv.size(); ++i) {
        ++symbols[iv[i]];
    }

    std::mt19937_64 rng;
    std::uniform_int_distribution<uint64_t> distribution(0, iv.size()-1);
    auto dice = bind(distribution, rng);
    for (uint64_t l = 0; l<cs.size(); ++l) {
        is[l] = dice();
        is[l] = is[l]%symbols[cs[l]]+1;
    }
}

template<class t_wt>
struct wt_trait {
    static uint64_t test_access(const t_wt& wt, const vector<size_type>& is, uint64_t mask, uint64_t times=100000000) {
        return ::test_access(wt, is, mask, times);
    }
    static uint64_t test_inverse_select(const t_wt& wt, const vector<size_type>& is, uint64_t mask, uint64_t times=100000000) {
        return ::test_inverse_select(wt, is, mask, times);
    }
};

template<class t_rac, class t_bitvector, class t_select, class t_select_zero>
struct wt_trait<wt_gmr_rs<t_rac, t_bitvector, t_select, t_select_zero>> {
    static uint64_t test_access(const wt_gmr_rs<t_rac, t_bitvector, t_select, t_select_zero>&, const vector<size_type>&, uint64_t, uint64_t) {
        return 0;
    }
    static uint64_t test_inverse_select(const wt_gmr_rs<t_rac, t_bitvector, t_select, t_select_zero>&, const vector<size_type>&, uint64_t, uint64_t) {
        return 0;
    }
};

// argv[1] = test case path  argv[2] = test case type  argv[3] test case name argv[4] wavelet tree id
int main(int argc, char* argv[])
{
    if (argc < 4) {
        cout << "Usage: file num_bytes testcase_name wt_id" << endl;
        return 1;
    }
    uint8_t type = argv[2][0]=='d' ? 'd' : argv[2][0]-'0';
    const uint64_t reps = 100000;
    uint64_t log_s = 20;
    uint64_t mask = (1<<log_s)-1;
    uint64_t check = 0;
    uint64_t size = 1<<log_s;

    // create values
    size_type k = 0;
    vector<value_type> cs(size);
    vector<size_type> is(size);
    vector<size_type> is2(size);
    vector<size_type> js(size);

    {
        int_vector<> iv;
        load_vector_from_file(iv, argv[1], type);
        random_cs(iv, cs);
        random_is_js(iv, is, js);
        prepare_for_select(iv, cs, is2);
    }
    // construct
    memory_monitor::start();
    auto start = timer::now();
    WT_TYPE wt;
    construct(wt, argv[1], type);
    auto stop = timer::now();
    memory_monitor::stop();
    cout << "# constructs_time = " << duration_cast<milliseconds>(stop-start).count()/(double)1000 << endl;
    cout << "# constructs_space = " << memory_monitor::peak() << endl;

    // size
    cout << "# wt_size = " << size_in_bytes(wt) << endl;

    // print structure
    // ofstream out("wt_"+string(argv[4])+"_"+string(argv[3])+".html");
    // write_structure<HTML_FORMAT>(wt, out);

    // access
    start = timer::now();
    check = wt_trait<WT_TYPE>::test_access(wt, is, mask, reps);
    stop = timer::now();
    cout << "# access_time = " << duration_cast<microseconds>(stop-start).count()/(double)reps << endl;
    cout << "# access_check = " << check << endl;

    // rank
    start = timer::now();
    check = test_rank(wt, is, cs, mask, reps);
    stop = timer::now();
    cout << "# rank_time = " << duration_cast<microseconds>(stop-start).count()/(double)reps << endl;
    cout << "# rank_check = " << check << endl;

    // inverse_select
    start = timer::now();
    check = wt_trait<WT_TYPE>::test_inverse_select(wt, is, mask, reps);
    stop = timer::now();
    cout << "# inverse_select_time = " << duration_cast<microseconds>(stop-start).count()/(double)reps << endl;
    cout << "# inverse_select_check = " << check << endl;

    // interval_symbols
    const uint64_t reps_interval_symbols = wt.sigma < 10000 ? reps : reps/100;
    start = timer::now();
    check = test_interval_symbols<WT_TYPE>(wt, is, js, k, mask, reps_interval_symbols);
    stop = timer::now();
    cout << "# interval_symbols_time = " << duration_cast<microseconds>(stop-start).count()/(double)reps_interval_symbols << endl;
    cout << "# interval_symbols_check = " << check << endl;

    // lex_count
    start = timer::now();
    check = test_lex_count<WT_TYPE>(wt, is, js, cs, mask, reps);
    stop = timer::now();
    cout << "# lex_count_time = " << duration_cast<microseconds>(stop-start).count()/(double)reps << endl;
    cout << "# lex_count_check = " << check << endl;

    // lex_smaller_count
    start = timer::now();
    check = test_lex_smaller_count<WT_TYPE>(wt, is, cs, mask, reps);
    stop = timer::now();
    cout << "# lex_smaller_count_time = " << duration_cast<microseconds>(stop-start).count()/(double)reps << endl;
    cout << "# lex_smaller_count_check = " << check << endl;

    // select
    start = timer::now();
    check = test_select(wt, is2, cs, mask, reps);
    stop = timer::now();
    cout << "# select_time = " << duration_cast<microseconds>(stop-start).count()/(double)reps << endl;
    cout << "# select_check = " << check << endl;

    return 0;
}
