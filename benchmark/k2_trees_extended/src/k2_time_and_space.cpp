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
    if (argc < 4) {
        cout << "Usage: input_graph queries k2tree_output" << endl;
        return 1;
    }


    K2_TYPE k2tree;
    k2_tree_ns::construction_algorithm construction = k2_tree_ns::ZORDER_SORT; //should be determined by type automatically
    bool use_shortcut = false;
    {
        mem_monitor mem_monitor1("bench_script_mem");
        auto start = timer::now();
        k2tree.load_from_ladrabin(argv[1], construction, 0, ram_file_name("asd"));
        auto stop = timer::now();
        auto status = mem_monitor1.get_current_stats();
        cout << "# constructs_time = " << duration_cast<milliseconds>(stop-start).count()/(double)1000 << endl;
        cout << "# constructs_space = " << status.VmHWM << endl;
        cout << "# constructs_space (VMEM) = " << status.VmPeak << endl;
    }

    // size
    cout << "# k2_size = " << size_in_bytes(k2tree) << endl;
    std::string output_file_name = argv[3];
    store_to_file(output_file_name, k2tree);


    access_times times_uncomp = perform_speed_test(argv[2], k2tree);


    cout << "# adj_time = " << times_uncomp.check_time << endl;
    cout << "# adj_check = " << times_uncomp.check_present << endl;

    cout << "# neighbors_time = " << times_uncomp.direct_time << endl;
    cout << "# neighbors_check = " << times_uncomp.direct_recovered << endl;

    cout << "# reverse_neighbors_time = " << times_uncomp.inverse_time << endl;
    cout << "# reverse_neighbors_check = " << times_uncomp.inverse_recovered << endl;

    {
        mem_monitor mem_monitor1("comp_mem");
        auto start = timer::now();
        k2tree.compress_leaves(DAC);
        auto stop = timer::now();
        auto status = mem_monitor1.get_current_stats();
        cout << "# compression_time = " << duration_cast<milliseconds>(stop-start).count()/(double)1000 << endl;
        cout << "# compression_space = " << status.VmHWM << endl;
        cout << "# compression_space (VMEM) = " << status.VmPeak << endl;
    }

    access_times times_comp = perform_speed_test(argv[2], k2tree);


    cout << "# adj_time_comp = " << times_comp.check_time << endl;
    cout << "# adj_check_comp = " << times_comp.check_present << endl;

    cout << "# neighbors_time_comp = " << times_comp.direct_time << endl;
    cout << "# neighbors_check_comp = " << times_comp.direct_recovered << endl;

    cout << "# reverse_neighbors_time_comp = " << times_comp.inverse_time << endl;
    cout << "# reverse_neighbors_check_comp = " << times_comp.inverse_recovered << endl;

    output_file_name.append("compressed");
    store_to_file(output_file_name, k2tree);

    return 0;
}
