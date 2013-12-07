#include <sdsl/suffix_trees.hpp>
#include <iostream>
#include <random>

using namespace sdsl;
using namespace std;

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

template<class t_csa>
void run(const char* file)
{
    mt19937_64 rng;
    rng.seed(424242);
    string tmp_csa = "tmp_csa.sdsl";
    {
        t_csa csa;
        construct(csa, file, 1);
        cout<<"csa_space_in_MB = "<<size_in_mega_bytes(csa)<<std::endl;
        cout<<"csa_sample_space_in_MB = "<<size_in_mega_bytes(csa.sa_sample)+size_in_mega_bytes(csa.isa_sample)<<std::endl;
        store_to_file(csa, tmp_csa);
    }
    t_csa csa;
    load_from_file(csa, tmp_csa);
    uniform_int_distribution<uint64_t> distribution(0, csa.size()-1);
    auto dice = bind(distribution, rng);

    uint64_t check=0;
    auto start = timer::now();
    size_t iterations = 1000000;
    for (size_t i=0; i<iterations; ++i) {
        check += csa.isa[dice()];
    }
    auto stop = timer::now();
    cout<<"check = "<<check<<endl;
    cout<<"time_in_us_per_isa_access  = "<< duration_cast<microseconds>(stop-start).count()/iterations << endl;

    check=0;
    start = timer::now();
    iterations = 1000000;
    for (size_t i=0; i<iterations; ++i) {
        check += csa[dice()];
    }
    stop = timer::now();
    cout<<"check = "<<check<<endl;
    cout<<"time_in_us_per_sa_access  = "<< duration_cast<microseconds>(stop-start).count()/iterations << endl;

    ofstream out("isa_access_"+to_string(t_csa::sa_sample_dens)+".html");
    write_structure<HTML_FORMAT>(csa, out);
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " file" << endl;
        cout << " Creates two CSAs for a byte file and times the ISA operation." << endl;
        return 1;
    }
    using t_wt = wt_huff<rrr_vector<63>>;
    const uint32_t S_SA =16;
    const uint32_t S_ISA=16;
    using t_sa_sample = fuzzy_sa_sampling<sd_vector<>>;
    run< csa_wt<t_wt,S_SA,S_ISA,t_sa_sample,fuzzy_isa_sampling_support<>> >(argv[1]);
    run< csa_wt<t_wt,2*S_SA,2*S_ISA,t_sa_sample,fuzzy_isa_sampling_support<>> >(argv[1]);
    run< csa_wt<t_wt,4*S_SA,4*S_ISA,t_sa_sample,fuzzy_isa_sampling_support<>> >(argv[1]);
    /*
        run< csa_wt<t_wt,S_SA,S_ISA,t_sa_sample,text_order_isa_sampling_support<inv_perm_support<8>>> >(argv[1]);
        run< csa_wt<t_wt,S_SA,S_ISA,t_sa_sample,text_order_isa_sampling_support<inv_perm_support<16>>> >(argv[1]);
        run< csa_wt<t_wt,S_SA,S_ISA,t_sa_sample,text_order_isa_sampling_support<inv_perm_support<32>>> >(argv[1]);
        run< csa_wt<t_wt,S_SA,S_ISA,t_sa_sample,text_order_isa_sampling_support<inv_perm_support<64>>> >(argv[1]);
    */
}
