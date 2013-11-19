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

//test access
template<class t_wt>
uint64_t test_access(const t_wt& wt, const vector<size_type>& is,uint64_t mask, uint64_t times=100000000)
{
    uint64_t cnt=0;
    for (uint64_t i=0; i<times; ++i) {
        cnt += wt[is[i&mask]];
    }
    return cnt;
}

//test rank
template<class t_wt>
uint64_t test_rank(const t_wt& wt, const vector<size_type>& is, const vector<value_type>& cs,uint64_t mask, uint64_t times=100000000)
{
    uint64_t cnt=0;
    for (uint64_t i=0; i<times; ++i) {
        cnt += wt.rank(is[i&mask],cs[i&mask]);
    }
    return cnt;
}

//test inverse_select
template<class t_wt>
uint64_t test_inverse_select(const t_wt& wt, const vector<size_type>& is, uint64_t mask, uint64_t times=100000000)
{
    uint64_t cnt=0;
    for (uint64_t i=0; i<times; ++i) {
        cnt += wt.inverse_select(is[i&mask]).first;
    }
    return cnt;
}

//test interval_symbols
template<class t_wt>
uint64_t test_interval_symbols(const t_wt& wt, const vector<size_type>& is, const vector<size_type>& js, size_type& k,vector<value_type>& tmp, vector<size_type>& tmp2, uint64_t mask, uint64_t times=100000000)
{
    uint64_t cnt=0;
    for (uint64_t i=0; i<times; ++i) {
        wt.interval_symbols(is[i&mask],js[i&mask],k,tmp,tmp2,tmp2);
        cnt += k;
    }
    return cnt;
}

//test lex_count
template<class t_wt>
uint64_t test_lex_count(const t_wt& wt, const vector<size_type>& is, const vector<size_type>& js, const vector<value_type>& cs, uint64_t mask, uint64_t times=100000000)
{
    uint64_t cnt=0;
    for (uint64_t i=0; i<times; ++i) {
        cnt += get<0>(wt.lex_count(is[i&mask],js[i&mask],cs[i&mask]));
    }
    return cnt;
}

//test lex_smaller_count
template<class t_wt>
uint64_t test_lex_smaller_count(const t_wt& wt, const vector<size_type>& is, const vector<value_type>& cs, uint64_t mask, uint64_t times=100000000)
{
    uint64_t cnt=0;
    for (uint64_t i=0; i<times; ++i) {
        cnt += get<0>(wt.lex_smaller_count(is[i&mask],cs[i&mask]));
    }
    return cnt;
}

//test select
template<class t_wt>
uint64_t test_select(const t_wt& wt,const vector<size_type>& is, const vector<value_type>& cs, uint64_t mask, uint64_t times=100000000)
{
    uint64_t cnt=0;
    for (uint64_t i=0; i<times; ++i) {
        cnt += wt.select(is[i&mask],cs[i&mask]);
    }
    return cnt;
}

//generate benchmark input
template<class t_wt>
void random_cs(const t_wt& wt,vector<value_type>& cs)
{

    std::mt19937_64 rng;
    std::uniform_int_distribution<uint64_t> distribution(0, wt.size()-1);
    auto dice = bind(distribution, rng);
    for (uint64_t l = 0; l<cs.size(); ++l) {
        cs[l]=wt[dice()];
    }
}

template<class t_wt>
void random_is_js(const t_wt& wt,vector<size_type>& is,vector<size_type>& js)
{

    std::mt19937_64 rng;
    std::uniform_int_distribution<uint64_t> distribution(0, wt.size()-1);
    auto dice = bind(distribution, rng);
    for (uint64_t l = 0; l<is.size(); ++l) {
        is[l]=dice();
        js[l]=min(is[l]+dice(),wt.size());
    }
}

template<class t_wt>
void prepare_for_select(const t_wt& wt,vector<value_type>& cs,vector<size_type>& is)
{

    std::mt19937_64 rng;
    std::uniform_int_distribution<uint64_t> distribution(0, wt.size()-1);
    auto dice = bind(distribution, rng);
    for (uint64_t l = 0; l<cs.size(); ++l) {
        is[l]=dice();
        is[l] = is[l]%(wt.rank(wt.size(),cs[l]))+1;
    }
}

template<class t_wt, bool lex_ordered = t_wt::lex_ordered>
struct wt_trait;

template<class t_wt>
struct wt_trait<t_wt, false> {
    static uint64_t test_interval_symbols(const t_wt& wt, const vector<size_type>& is, const vector<size_type>& js, size_type& k,vector<value_type>& tmp, vector<size_type>& tmp2, uint64_t mask, uint64_t times=100000000) {
        return ::test_interval_symbols(wt,is,js,k,tmp,tmp2,mask,times);
    }
    static uint64_t test_lex_count(const t_wt& wt, const vector<size_type>& is, const vector<size_type>& js, const vector<value_type>& cs, uint64_t mask, uint64_t times=100000000) {
        return 0;
    }
    static uint64_t test_lex_smaller_count(const t_wt& wt, const vector<size_type>& is, const vector<value_type>& cs, uint64_t mask, uint64_t times=100000000) {
        return 0;
    }
};

template<class t_wt>
struct wt_trait<t_wt, true> {
    static uint64_t test_interval_symbols(const t_wt& wt, const vector<size_type>& is, const vector<size_type>& js, size_type& k,vector<value_type>& tmp, vector<size_type>& tmp2, uint64_t mask, uint64_t times=100000000) {
        return ::test_interval_symbols(wt,is,js,k,tmp,tmp2,mask,times);
    }
    static uint64_t test_lex_count(const t_wt& wt, const vector<size_type>& is, const vector<size_type>& js, const vector<value_type>& cs, uint64_t mask, uint64_t times=100000000) {
        return ::test_lex_count(wt,is,js,cs,mask,times);
    }
    static uint64_t test_lex_smaller_count(const t_wt& wt, const vector<size_type>& is, const vector<value_type>& cs, uint64_t mask, uint64_t times=100000000) {
        return ::test_lex_smaller_count(wt,is,cs,mask,times);
    }
};

template<class t_bitvector, class t_rank, class t_select, class t_wt>
struct wt_trait<wt_rlmn<t_bitvector, t_rank, t_select, t_wt>, false> {
    static uint64_t test_interval_symbols(const wt_rlmn<t_bitvector, t_rank, t_select, t_wt>& wt, const vector<size_type>& is, const vector<size_type>& js, size_type& k,vector<value_type>& tmp, vector<size_type>& tmp2, uint64_t mask, uint64_t times=100000000) {
        return 0;
    }
    static uint64_t test_lex_count(const wt_rlmn<t_bitvector, t_rank, t_select, t_wt>& wt, const vector<size_type>& is, const vector<size_type>& js, const vector<value_type>& cs, uint64_t mask, uint64_t times=100000000) {
        return 0;
    }
    static uint64_t test_lex_smaller_count(const wt_rlmn<t_bitvector, t_rank, t_select, t_wt>& wt, const vector<size_type>& is, const vector<value_type>& cs, uint64_t mask, uint64_t times=100000000) {
        return 0;
    }
};

// argv[1] = test case path  argv[2] = test case type
int main(int argc, char* argv[])
{
    uint8_t type = argv[2][0]=='d' ? 'd' : argv[2][0]-'0';

    WT_TYPE wt;
    //construct
    auto start = timer::now();
    construct(wt,argv[1],type);
    auto stop = timer::now();
    cout << "# constructs_time = " << duration_cast<seconds>(stop-start).count() << endl;

    //size
    cout << "# wt_size = " << size_in_bytes(wt) << endl;

    const uint64_t reps = 100000;
    uint64_t log_s = 20;
    uint64_t mask = (1<<log_s)-1;
    uint64_t check = 0;
    uint64_t size = 1<<log_s;

    //create values
    size_type k =0;
    vector<value_type> cs(size);
    vector<size_type> is(size);
    vector<size_type> js(size);
    vector<value_type> tmp(wt.sigma);
    vector<size_type> tmp2(wt.sigma);
    random_cs(wt,cs);
    random_is_js(wt,is,js);

    //access
    start = timer::now();
    check = test_access(wt,is,mask,reps);
    stop = timer::now();
    cout << "# access_time = " << duration_cast<microseconds>(stop-start).count()/(double)reps << endl;
    cout << "# access_check = " << check << endl;

    //rank
    start = timer::now();
    check = test_rank(wt,is,cs,mask,reps);
    stop = timer::now();
    cout << "# rank_time = " << duration_cast<microseconds>(stop-start).count()/(double)reps << endl;
    cout << "# rank_check = " << check << endl;

    //inverse_select
    start = timer::now();
    check = test_inverse_select(wt,is,mask,reps);
    stop = timer::now();
    cout << "# inverse_select_time = " << duration_cast<microseconds>(stop-start).count()/(double)reps << endl;
    cout << "# inverse_select_check = " << check << endl;

    //interval_symbols
    start = timer::now();
    check = wt_trait<WT_TYPE>::test_interval_symbols(wt,is,js,k,tmp,tmp2,mask,reps);
    stop = timer::now();
    cout << "# interval_symbols_time = " << duration_cast<microseconds>(stop-start).count()/(double)reps << endl;
    cout << "# interval_symbols_check = " << check << endl;

    //lex_count
    start = timer::now();
    check = wt_trait<WT_TYPE>::test_lex_count(wt,is,js,cs,mask,reps);
    stop = timer::now();
    cout << "# lex_count_time = " << duration_cast<microseconds>(stop-start).count()/(double)reps << endl;
    cout << "# lex_count_check = " << check << endl;

    //lex_smaller_count
    start = timer::now();
    check = wt_trait<WT_TYPE>::test_lex_smaller_count(wt,is,cs,mask,reps);
    stop = timer::now();
    cout << "# lex_smaller_count_time = " << duration_cast<microseconds>(stop-start).count()/(double)reps << endl;
    cout << "# lex_smaller_count_check = " << check << endl;

    prepare_for_select(wt,cs,is);

    //select
    start = timer::now();
    check = test_select(wt,is,cs,mask,reps);
    stop = timer::now();
    cout << "# select_time = " << duration_cast<microseconds>(stop-start).count()/(double)reps << endl;
    cout << "# select_check = " << check << endl;
}
