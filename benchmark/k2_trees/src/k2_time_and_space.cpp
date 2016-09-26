#include<fstream>
#include<string>
#include<sstream>

#include<sdsl/bit_vectors.hpp>
#include<sdsl/k2_tree.hpp>

using namespace std;
using namespace sdsl;
using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

typedef K2_TYPE::idx_type idx_type;
typedef K2_TYPE::size_type size_type;


template<class t_kt>
size_type test_adj(const t_kt &tree, idx_type node, size_type neighbor,
				   uint64_t times)
{
	size_type cnt = 0;
	for(uint64_t i = 0; i < times; i++)
		if(tree.adj(node, neighbor))
				cnt++;
	return cnt;
}

template<class t_kt>
size_type test_neighbors(const t_kt &tree, idx_type node, uint64_t times)
{
	size_type cnt = 0;
	for(uint64_t i = 0; i < times; i++)
		cnt += tree.neigh(node).size();
	return cnt;
}

template<class t_kt>
size_type test_reverse_neighbors(const t_kt &tree, idx_type node, uint64_t times)
{
	size_type cnt = 0;
	for(uint64_t i = 0; i < times; i++)
		cnt += tree.reverse_neigh(node).size();
	return cnt;
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        cout << "Usage: file" << endl;
        return 1;
    }

    const uint64_t reps = 100000;

    // construct
    memory_monitor::start();
	std::ifstream is(argv[1]);
    auto start = timer::now();
    K2_TYPE k2;
    k2.load(is);
    auto stop = timer::now();
    memory_monitor::stop();
    cout << "# constructs_time = " << duration_cast<milliseconds>(stop-start).count()/(double)1000 << endl;
    cout << "# constructs_space = " << memory_monitor::peak() << endl;
    // size
    cout << "# k2_size = " << size_in_bytes(k2) << endl;
    is.close();

	// adj
	start = timer::now();
	auto check = test_adj<K2_TYPE>(k2, 7, 5, reps);
    stop = timer::now();
    cout << "# adj_time = " << duration_cast<microseconds>(stop-start).count()/(double)reps << endl;
    cout << "# adj_check = " << check << endl;

    // neighbors
	start = timer::now();
	check = test_neighbors<K2_TYPE>(k2, 7, reps);
    stop = timer::now();
    cout << "# neighbors_time = " << duration_cast<microseconds>(stop-start).count()/(double)reps << endl;
    cout << "# neighbors_check = " << check << endl;

    start = timer::now();
	check = test_reverse_neighbors<K2_TYPE>(k2, 10, reps);
    stop = timer::now();
    cout << "# reverse_neighbors_time = " << duration_cast<microseconds>(stop-start).count()/(double)reps << endl;
    cout << "# reverse_neighbors_check = " << check << endl;

    return 0;
}
