#include<fstream>
#include<string>
#include<sstream>

#include<sdsl/bit_vectors.hpp>
#include<sdsl/k2_tree_algorithm.hpp>
#include <sdsl/mem_monitor.hpp>
#include <sdsl/k2_tree_utility.hpp>
#include <sdsl/k2_tree_compressor.hpp>

using namespace std;
using namespace sdsl;
using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

int main(int argc, char* argv[])
{
    if (argc < 5) {
        cout << "Usage: input_graph queries k2tree_output split_size" << endl;
        return 1;
    }

    std::string query_file = argv[2];
    K2_TYPE k2tree;
    k2_tree_ns::construction_algorithm construction = k2_tree_ns::ZORDER_SORT; //should be determined by type automatically
    bool use_shortcut = false;
    uint64_t construction_time;
    uint64_t peak_RSS;
    uint64_t peak_VMEM;
    {
        mem_monitor mem_monitor1("bench_script_mem");
        auto start = timer::now();
        cout << "matrix_shift = " << atoi(argv[4]) << endl;
        k2tree.load_from_ladrabin(argv[1], construction, "", 0, atoi(argv[4]));
        auto stop = timer::now();
        auto status = mem_monitor1.get_current_stats();
        construction_time = duration_cast<seconds>(stop-start).count();
        cout << "# constructs_time = " << duration_cast<seconds>(stop-start).count() << endl;
        cout << "# constructs_space = " << status.VmHWM << endl;
        cout << "# constructs_space_vmem = " << status.VmPeak << endl;
        peak_RSS = status.VmHWM;
        peak_VMEM = status.VmPeak;

        cout << "# construct_morton_duration = " << morton_number_duration << endl;
        cout << "# construct_bv_complete_duration = " << construct_bv_complete_duration  << endl;
        cout << "# construct_sort_duration = " << sort_duration << endl;
        cout << "# construct_duration = " << construct_duration << endl;
        cout << "# buildvec_duration = " << build_vec_duration << endl;
        cout << "# subtree_constructor_duration = " << constructor_duration << endl;
        cout << "# constructor_call_duration = " << constructor_call_duration << endl;
    }

    // size
    cout << "# uncompressed_size = " << size_in_bytes(k2tree) << endl;
    std::string output_file_name = argv[3];
    store_to_file(output_file_name, k2tree);

    access_times times_uncomp = perform_speed_test(query_file, k2tree, query_file+".single");


    cout << "# adj_time = " << times_uncomp.check_time << endl;
    cout << "# adj_check = " << times_uncomp.check_present << endl;

    cout << "# neighbors_time = " << times_uncomp.direct_time << endl;
    cout << "# neighbors_check = " << times_uncomp.direct_recovered << endl;

    cout << "# reverse_neighbors_time = " << times_uncomp.inverse_time << endl;
    cout << "# reverse_neighbors_check = " << times_uncomp.inverse_recovered << endl;

    std::cout << " CSV = " << argv[1] << "\t" << k2tree.get_type_string() <<  "\t" << construction_time << "\t" << peak_RSS << "\t" << peak_VMEM;
    if (use_shortcut){
        //Construction Time	Compressed Size (Byte)	Bpe	Direct Short (ns)	Direct (ns)	Inverse Short (ns)	Inverse (ns)	Check S (ns)	Check (ns)
        std::cout << "\t" <<  times_uncomp.direct_short_time << "\t" << times_uncomp.direct_time << "\t" << times_uncomp.inverse_short_time << "\t" << times_uncomp.inverse_time << "\t" << times_uncomp.check_short_time << "\t" << times_uncomp.check_time << std::endl;
    } else {
        //Construction Time	Compressed Size (Byte)	Bpe	Direct (ns)	Inverse (ns)	Check (ns)
        std::cout << "\t" << times_uncomp.direct_time << "\t" << times_uncomp.inverse_time << "\t" << times_uncomp.check_time << std::endl;
    }

    uint64_t peak_RSS_comp;
    uint64_t peak_VMEM_comp;
    uint64_t construction_time_comp;
    {
        mem_monitor mem_monitor1("comp_mem");
        auto start = timer::now();
        k2tree.compress_leaves(DAC);
        auto stop = timer::now();
        auto status = mem_monitor1.get_current_stats();
        cout << "# compression_time = " << duration_cast<seconds>(stop-start).count() << endl;
        cout << "# compression_space = " << status.VmHWM << endl;
        cout << "# compression_space_vmem = " << status.VmPeak << endl;
        peak_RSS_comp = mem_monitor1.max_seen_rss;
        peak_VMEM_comp = mem_monitor1.max_seen_vmem;
        construction_time_comp = duration_cast<milliseconds>(stop - start).count();
    }

    cout << "# compressed_size = " << size_in_bytes(k2tree) << endl;

    access_times times_comp = perform_speed_test(argv[2], k2tree, query_file+".single");


    cout << "# adj_time_comp = " << times_comp.check_time << endl;
    cout << "# adj_check_comp = " << times_comp.check_present << endl;

    cout << "# neighbors_time_comp = " << times_comp.direct_time << endl;
    cout << "# neighbors_check_comp = " << times_comp.direct_recovered << endl;

    cout << "# reverse_neighbors_time_comp = " << times_comp.inverse_time << endl;
    cout << "# reverse_neighbors_check_comp = " << times_comp.inverse_recovered << endl;

    cout << "# comp_word_it = " << word_iteration << endl;
    cout << "# frequency_encoding = " << frequency_encoding << endl;
    cout << "# dac_compression = " << dac_compression << endl;

    output_file_name.append("compressed");
    store_to_file(output_file_name, k2tree);

    std::cout << " CSV = " << argv[1] << "\t" << k2tree.get_type_string() <<  "\t" << construction_time_comp << "\t" << peak_RSS_comp << "\t" << peak_VMEM_comp;
    if (use_shortcut){
        //Construction Time	Compressed Size (Byte)	Bpe	Direct Short (ns)	Direct (ns)	Inverse Short (ns)	Inverse (ns)	Check S (ns)	Check (ns)
        std::cout << "\t" << times_comp.direct_short_time <<"\t"<< times_comp.direct_time <<"\t"<< times_comp.inverse_short_time <<"\t"<< times_comp.inverse_time <<"\t"<< times_comp.check_short_time <<"\t"<< times_comp.check_time << std::endl;
    } else {
        //Construction Time	Compressed Size (Byte)	Bpe	Direct (ns)	Inverse (ns)	Check (ns)
        std::cout << "\t" << times_comp.direct_time << "\t" << times_comp.inverse_time << "\t" << times_comp.check_time << std::endl;
    }

    return 0;
}
