#include <chrono>
#include "sdsl/k2_tree_algorithm.hpp"
#include "gtest/gtest.h"

namespace {

    using namespace sdsl;
    using namespace std;
    using namespace std::chrono;
    using timer = std::chrono::high_resolution_clock;

    typedef int_vector<>::size_type size_type;

    string ladrabin_file;
    string temp_file;
    string query_file;
    uint queryCount;
    uint *queries;
    uint construction_time = 0;
    uint compressed_size = 0;

    template<class T>
    class k2_performance_test : public ::testing::Test {
    };


// Writes primitive-typed variable t to stream out
    template<class T>
    size_t _write_member(const T &t, std::ostream &out) {
        out.write((char *) &t, sizeof(t));
        size_t written_bytes = sizeof(t);
        return written_bytes;
    }

// Writes array-typed variable t of length length to stream out
    template<class T>
    size_t _write_member(const T *t, size_t length, std::ostream &out) {
        out.write((char *) t, length * sizeof(T));
        size_t written_bytes = length * sizeof(T);
        return written_bytes;
    }

// Reads primitive-typed variable t from stream om
    template<class T>
    void _read_member(T &t, std::istream &in) {
        in.read((char *) &t, sizeof(t));
    }

// Reads array-typed variable t of length length from stream in
    template<class T>
    void _read_member(T **t, size_t length, std::istream &in) {
        *t = (T *) malloc(sizeof(T) * length);
        in.read((char *) *t, length * sizeof(T));
    }

    int get_file_size(std::string filename) // path to file
    {
        FILE *p_file = NULL;
        p_file = fopen(filename.c_str(), "rb");
        fseek(p_file, 0, SEEK_END);
        int size = ftell(p_file);
        fclose(p_file);
        return size;
    }

    using testing::Types;

    typedef k2_tree_hybrid<2, 2, 2, 2, bit_vector, bit_vector, true> hybrid_k2_2222_b_comp;
    typedef k2_tree_hybrid<4, 5, 2, 4, bit_vector, bit_vector, true> hybrid_k2_4524_b_comp;
    //typedef k2_tree_hybrid<8, 5, 2, 4, bit_vector, bit_vector, true> hybrid_k2_8524_b_comp;
    typedef k2_tree_hybrid<8, 5, 2, 8, bit_vector, bit_vector, true> hybrid_k2_8528_b_comp;
    typedef k2_tree_hybrid<8, 3, 2, 8, bit_vector, bit_vector, true> hybrid_k2_8328_b_comp;
    //typedef k2_tree_hybrid<8, 3, 2, 4, bit_vector, bit_vector, true> hybrid_k2_8324_b_comp;
    typedef k2_tree_hybrid<4, 5, 2, 8, bit_vector, bit_vector, true> hybrid_k2_4528_b_comp;
    typedef k2_tree_hybrid<4, 6, 2, 4, bit_vector, bit_vector, true> hybrid_k2_4624_b_comp;
    typedef k2_tree_hybrid<2, 5, 2, 8, bit_vector, bit_vector, true> hybrid_k2_2528_b_comp;
    typedef k2_tree_hybrid<16, 5, 2, 16, bit_vector, bit_vector, true> hybrid_k2_165216_b_comp;
    //typedef k2_tree_hybrid<8, 5, 2, 8, bit_vector, bit_vector, true> hybrid_k2_8523_b_comp;
    typedef k2_tree_hybrid<3, 5, 2, 4, bit_vector, bit_vector, true> hybrid_k2_3524_b_comp;
    //typedef k2_tree_hybrid<2, 2, 2, 2, bit_vector, bit_vector, false> hybrid_k2_2222_b;
    typedef k2_tree_hybrid<4, 5, 2, 4, bit_vector, bit_vector, false> hybrid_k2_4524_b;
    typedef k2_tree_hybrid<4, 6, 2, 4, bit_vector, bit_vector, false> hybrid_k2_4624_b;
    typedef k2_tree_hybrid<2, 5, 2, 8, bit_vector, bit_vector, false> hybrid_k2_2528_b;
    typedef k2_tree_hybrid<16, 5, 2, 16, bit_vector, bit_vector, false> hybrid_k2_165216_b;
    typedef k2_tree_hybrid<8, 5, 2, 8, bit_vector, bit_vector, false> hybrid_k2_8523_b;
    //typedef k2_tree_hybrid<3, 5, 2, 4, bit_vector, bit_vector, false> hybrid_k2_3524_b;
    //typedef k2_tree_hybrid<8, 5, 2, 4, bit_vector, bit_vector, false> hybrid_k2_8524_b;
    typedef k2_tree_hybrid<8, 5, 2, 8, bit_vector, bit_vector, false> hybrid_k2_8528_b;
    typedef k2_tree_hybrid<8, 3, 2, 8, bit_vector, bit_vector, false> hybrid_k2_8328_b;
    //typedef k2_tree_hybrid<8, 3, 2, 4, bit_vector, bit_vector, false> hybrid_k2_8324_b;
    typedef k2_tree_hybrid<4, 5, 2, 8, bit_vector, bit_vector, false> hybrid_k2_4528_b;


    typedef k2_tree<2, bit_vector> k2;
    typedef k2_tree<2, rrr_vector<63>> k2rrr;
    typedef k2_tree<3, bit_vector> k3;
    typedef k2_tree<4, bit_vector> k4;
    typedef k2_tree<6, bit_vector> k6;
    typedef k2_tree<8, bit_vector> k8;
    typedef k2_tree<16, bit_vector> k16;

    typedef k2_tree<2, bit_vector, bit_vector, true> k2_comp;
    typedef k2_tree<3, bit_vector, bit_vector, true> k3_comp;
    typedef k2_tree<4, bit_vector, bit_vector, true> k4_comp;
    typedef k2_tree<6, bit_vector, bit_vector, true> k6_comp;
    typedef k2_tree<8, bit_vector, bit_vector, true> k8_comp;
    typedef k2_tree<16, bit_vector, bit_vector, true> k16_comp;

    typedef k2_tree<2, bit_vector, bit_vector, true, 2> k2_comp_2;
    typedef k2_tree<4, bit_vector, bit_vector, true, 2> k4_comp_2;
    typedef k2_tree<8, bit_vector, bit_vector, true, 2> k8_comp_2;
    typedef k2_tree<16, bit_vector, bit_vector, true, 2> k16_comp_2;

    typedef k2_tree<2, bit_vector, bit_vector, true, 4> k2_comp_4;
    typedef k2_tree<4, bit_vector, bit_vector, true, 4> k4_comp_4;
    typedef k2_tree<8, bit_vector, bit_vector, true, 4> k8_comp_4;
    typedef k2_tree<16, bit_vector, bit_vector, true, 4> k16_comp_4;


    typedef k2_tree<8, bit_vector, bit_vector, true, 5> k8_comp_5;
    //typedef k2_tree<2, bit_vector, bit_vector, true, 8> k2_comp_8;
    //typedef k2_tree<4, bit_vector, bit_vector, true, 8> k4_comp_8;


    typedef k2_tree<2, bit_vector, bit_vector, false, 2> k2_2;
    typedef k2_tree<4, bit_vector, bit_vector, false, 2> k4_2;
    typedef k2_tree<8, bit_vector, bit_vector, false, 2> k8_2;
    typedef k2_tree<16, bit_vector, bit_vector, false, 2> k16_2;

    typedef k2_tree<2, bit_vector, bit_vector, false, 4> k2_4;
    //typedef k2_tree<4, bit_vector, bit_vector, false, 4> k4_4;
    //typedef k2_tree<8, bit_vector, bit_vector, false, 4> k8_4;
    //typedef k2_tree<16, bit_vector, bit_vector, false, 4> k16_4;


    //typedef k2_tree<8, bit_vector, bit_vector, false, 5> k8_5;
    typedef k2_tree<2, bit_vector, bit_vector, false, 8> k2_8;
    //typedef k2_tree<4, bit_vector, bit_vector, false, 8> k4_8;

    typedef Types<
            /*k2,
            k3,
            k4,
            k6,
            k8,
            k16,*/
/*            k2_comp,
            k3_comp,
            k4_comp,
            k6_comp,
            k8_comp,
            k16_comp,
            k2_comp_2,
            k4_comp_2,
            k8_comp_2,
            k16_comp_2,
            k2_comp_4,
            k4_comp_4,
            k8_comp_4,
            k16_comp_4,
            k8_comp_5,
            k2_tree_partitioned<4,k2,true>,
            k2_tree_partitioned<4,k3,true>,
            k2_tree_partitioned<4,k4,true>,
            k2_tree_partitioned<4,k6,true>,
            k2_tree_partitioned<4,k8,true>,
            k2_tree_partitioned<4,k16,true>,
            k2_tree_partitioned<4,k2_2,true>,
            k2_tree_partitioned<4,k4_2,true>,
            k2_tree_partitioned<4,k8_2,true>,
            k2_tree_partitioned<4,k16_2,true>,

            k2_tree_partitioned<8,k2,true>,
            k2_tree_partitioned<8,k3,true>,
            k2_tree_partitioned<8,k4,true>,
            k2_tree_partitioned<8,k6,true>,
            k2_tree_partitioned<8,k8,true>,
            k2_tree_partitioned<8,k16,true>,
            k2_tree_partitioned<8,k2_2,true>,
            k2_tree_partitioned<8,k4_2,true>,
            k2_tree_partitioned<8,k8_2,true>,
            k2_tree_partitioned<8,k16_2,true>,

            k2_tree_partitioned<16,k2,true>,
            k2_tree_partitioned<16,k3,true>,
            k2_tree_partitioned<16,k4,true>,
            k2_tree_partitioned<16,k6,true>,
            k2_tree_partitioned<16,k8,true>,
            k2_tree_partitioned<16,k16,true>,
            k2_tree_partitioned<16,k2_2,true>,
            k2_tree_partitioned<16,k4_2,true>,
            k2_tree_partitioned<16,k8_2,true>,
            k2_tree_partitioned<16,k16_2,true>,*/


            hybrid_k2_4524_b_comp,
            hybrid_k2_4528_b_comp,
            hybrid_k2_8528_b_comp,
            //hybrid_k2_8324_b_comp,
            //hybrid_k2_8328_b_comp,
            hybrid_k2_4624_b_comp,
            hybrid_k2_2528_b_comp,
            hybrid_k2_165216_b_comp,
            //hybrid_k2_8523_b_comp,
            hybrid_k2_3524_b_comp,

            k2_tree_partitioned<4,hybrid_k2_4524_b,true>,
            k2_tree_partitioned<4,hybrid_k2_4528_b,true>,
            k2_tree_partitioned<4,hybrid_k2_8528_b,true>,
            k2_tree_partitioned<4,hybrid_k2_8328_b,true>,
            k2_tree_partitioned<4,hybrid_k2_4624_b,true>,
            k2_tree_partitioned<4,hybrid_k2_2528_b,true>,
            k2_tree_partitioned<4,hybrid_k2_165216_b,true>,
            k2_tree_partitioned<4,hybrid_k2_8523_b,true>,

            k2_tree_partitioned<8,hybrid_k2_4524_b,true>,
            k2_tree_partitioned<8,hybrid_k2_4528_b,true>,
            k2_tree_partitioned<8,hybrid_k2_8528_b,true>,
            k2_tree_partitioned<8,hybrid_k2_8328_b,true>,
            k2_tree_partitioned<8,hybrid_k2_4624_b>,
            k2_tree_partitioned<8,hybrid_k2_2528_b>,
            k2_tree_partitioned<8,hybrid_k2_165216_b>,
            k2_tree_partitioned<8,hybrid_k2_8523_b>
    > Implementations;

    TYPED_TEST_CASE(k2_performance_test, Implementations);

    TYPED_TEST(k2_performance_test, CreateAndStoreTest) {
        auto start = timer::now();
        std::fstream fileStream(ladrabin_file, std::ios_base::in);
        if (fileStream.is_open()) {

            uint number_of_nodes;
            ulong number_of_edges;
            _read_member(number_of_nodes, fileStream);
            _read_member(number_of_edges, fileStream);

            uint nodes_read = 0;
            uint edges_read = 0;
            uint source_id;
            int target_id;

            std::vector<std::pair<uint, uint>> coords(number_of_edges);
            for (uint64_t i = 0; i < number_of_nodes + number_of_edges; i++) {
                _read_member(target_id, fileStream);
                if (target_id < 0) {
                    nodes_read++;
                } else {
                    source_id = nodes_read - 1;
                    coords.push_back(std::make_pair(source_id, target_id));
                    edges_read++;
                }
            }
            fileStream.close();

            TypeParam k2tree("", true, coords, number_of_nodes - 1);
            coords.clear();

            store_to_file(k2tree, temp_file);
            auto stop = timer::now();

            std::cout << "Construction time: " << duration_cast<milliseconds>(stop - start).count() << "\n";

            construction_time = duration_cast<milliseconds>(stop - start).count();
            std::cout << "Compressed size: " << get_file_size(temp_file) << std::endl;
            compressed_size = get_file_size(temp_file);
        } else {
            throw "Could not load file";
        }
    }


    TYPED_TEST(k2_performance_test, load_queries) {
        FILE *list_fp = fopen(query_file.c_str(), "r");
        auto tmp = fread(&queryCount, sizeof(uint), 1, list_fp);
        queries = (uint *) malloc(sizeof(uint) * queryCount);
        auto tmp2 = fread(queries, sizeof(uint), queryCount, list_fp);
        fclose(list_fp);
    }

    TYPED_TEST(k2_performance_test, direct_links_test) {
        TypeParam k2tree;
        load_from_file(k2tree, temp_file);

        bool use_shortut = true;
        auto direct_short_time = 0;
        auto direct_time = 0;
        auto inverse_short_time = 0;
        auto inverse_time = 0;
        auto check_short_time= 0;
        auto check_time= 0;

        try {
            std::vector<uint32_t> result;
            k2tree.direct_links_shortcut_2((uint) 0, result);
        } catch (std::runtime_error const &e) {
            use_shortut = false;
            std::cout << "Checking with shortcut deactivated, none present";
        }
        {
            uint i;
            ulong recovered = 0;
            std::vector<uint32_t> result;
            if (use_shortut) {
                std::cout << "Direct with access shortcut2" << std::endl;
                auto start = timer::now();
                for (i = 0; i < queryCount; i++) {
                    k2tree.direct_links_shortcut_2(queries[i], result);
                    recovered += result.size();
                }
                auto stop = timer::now();

                std::cout << "Recovered Nodes:" << recovered << "\n";
                std::cout << "Queries:" << queryCount << "\n";
                std::cout << "Total time(ns): " << duration_cast<nanoseconds>(stop - start).count() << "\n";
                std::cout << "Time per query: " << duration_cast<nanoseconds>(stop - start).count() / queryCount << "\n";
                std::cout << "Time per link: " << duration_cast<nanoseconds>(stop - start).count() / recovered << "\n";

                direct_short_time = duration_cast<nanoseconds>(stop - start).count() /recovered;
            }
        }

        {
            uint i;
            ulong recovered = 0;
            std::vector<uint32_t> result;
            std::cout << "Direct without access shortcut" << std::endl;
            auto start = timer::now();
            for (i = 0; i < queryCount; i++) {
                k2tree.direct_links2(queries[i], result);
                recovered += result.size();
            }
            auto stop = timer::now();

            std::cout << "Recovered Nodes:" << recovered << "\n";
            std::cout << "Queries:" << queryCount << "\n";
            std::cout << "Total time(ns): " << duration_cast<nanoseconds>(stop - start).count() << "\n";
            std::cout << "Time per query: " << duration_cast<nanoseconds>(stop - start).count() / queryCount << "\n";
            std::cout << "Time per link: " << duration_cast<nanoseconds>(stop - start).count() / recovered << "\n";

            direct_time = duration_cast<nanoseconds>(stop - start).count() /recovered;
        }

        {
            uint i;
            ulong recovered = 0;
            std::vector<uint32_t> result;
            if (use_shortut) {
                std::cout << "Inverse with access shortcut" << std::endl;
                auto start = timer::now();
                for (i = 0; i < queryCount; i++) {
                    k2tree.inverse_links_shortcut(queries[i], result);
                    recovered += result.size();
                }
                auto stop = timer::now();

                std::cout << "Recovered Nodes:" << recovered << "\n";
                std::cout << "Queries:" << queryCount << "\n";
                std::cout << "Total time(ns): " << duration_cast<nanoseconds>(stop - start).count() << "\n";
                std::cout << "Time per query: " << duration_cast<nanoseconds>(stop - start).count() / queryCount << "\n";
                std::cout << "Time per link: " << duration_cast<nanoseconds>(stop - start).count() / recovered << "\n";

                inverse_short_time = duration_cast<nanoseconds>(stop - start).count() / recovered;
            }
        }

        {
            uint i;
            ulong recovered = 0;
            std::vector<uint32_t> result;
            std::cout << "Inverse without access shortcut" << std::endl;
            auto start = timer::now();
            for (i = 0; i < queryCount; i++) {
                k2tree.inverse_links2(queries[i], result);
                recovered += result.size();
            }
            auto stop = timer::now();

            std::cout << "Recovered Nodes:" << recovered << "\n";
            std::cout << "Queries:" << queryCount << "\n";
            std::cout << "Total time(ns): " << duration_cast<nanoseconds>(stop - start).count() << "\n";
            std::cout << "Time per query: " << duration_cast<nanoseconds>(stop - start).count() / queryCount << "\n";
            std::cout << "Time per link: " << duration_cast<nanoseconds>(stop - start).count() / recovered << "\n";

            inverse_time = duration_cast<nanoseconds>(stop - start).count() / recovered;
        }

        srand(0);
        uint link_query_count = 1000000;
        std::vector<std::pair<uint, uint>> check_link_queries(link_query_count);
        for (uint i = 0; i < link_query_count; i++) {
            check_link_queries.push_back(
                    std::make_pair(rand() % (k2tree.m_max_element - 1), rand() % (k2tree.m_max_element - 1)));
        }

        {
            if (use_shortut) {
                std::cout << "Performing single link with shortcut" << std::endl;
                auto start = timer::now();
                for (auto pair: check_link_queries) {
                    k2tree.check_link_shortcut(pair);
                }
                auto stop = timer::now();

                std::cout << "Queries:" << link_query_count << "\n";
                std::cout << "Total time(ns): " << duration_cast<nanoseconds>(stop - start).count() << "\n";
                std::cout << "Time per query: " << duration_cast<nanoseconds>(stop - start).count() / link_query_count
                          << "\n";

                check_short_time = duration_cast<nanoseconds>(stop - start).count() / link_query_count;
            }
        }

        {
            std::cout << "Performing single link without shortcut" << std::endl;
            auto start = timer::now();
            for (auto pair: check_link_queries) {
                k2tree.check_link(pair);
            }
            auto stop = timer::now();

            std::cout << "Queries:" << link_query_count << "\n";
            std::cout << "Total time(ns): " << duration_cast<nanoseconds>(stop - start).count() << "\n";
            std::cout << "Time per query: " << duration_cast<nanoseconds>(stop - start).count() / link_query_count << "\n";

            check_time = duration_cast<nanoseconds>(stop - start).count()/ link_query_count;
        }


        if (use_shortut){
            //Construction Time	Compressed Size (Byte)	Bpe	Direct Short (ns)	Direct (ns)	Inverse Short (ns)	Inverse (ns)	Check S (ns)	Check (ns)
            std::cout << "Hereyougo:" << construction_time << "," << compressed_size << "," << direct_short_time <<","<< direct_time <<","<< inverse_short_time <<","<< inverse_time <<","<< check_short_time <<","<< check_time << std::endl;
        } else {
            //Construction Time	Compressed Size (Byte)	Bpe	Direct (ns)	Inverse (ns)	Check (ns)
            std::cout << "Hereyougo:" << construction_time << "," << compressed_size << "," << direct_time <<"," << inverse_time <<"," << check_time << std::endl;
        }
    }
}  // namespace

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    if (argc < 3) {
        // LCOV_EXCL_START
        cout << "Usage: " << argv[0] << " ladrabin_file temp_file query_file" << endl;
        cout << " (1) Generates a k2-tree out of ladrabin_file" << endl;
        cout << "     Result is stored in temp_file." << endl;
        cout << " (2) Performs performance tests." << endl;
        cout << " (3) Deletes temp_file." << endl;
        return 1;
        // LCOV_EXCL_STOP
    }
    ladrabin_file = argv[1];
    temp_file = argv[2];
    query_file = argv[3];

    return RUN_ALL_TESTS();
}
