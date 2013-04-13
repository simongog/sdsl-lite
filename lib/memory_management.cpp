#include "sdsl/memory_management.hpp"

#include <cstdlib> // for malloc and free
#include <sys/mman.h>

#ifdef MAP_HUGETLB
#define HUGE_LEN 1073741824
#define HUGE_PROTECTION (PROT_READ | PROT_WRITE)
#define HUGE_FLAGS (MAP_HUGETLB | MAP_ANONYMOUS | MAP_PRIVATE)
#endif

static int nifty_counter = 0;

std::map<uint64_t, sdsl::mm_item_base*> sdsl::mm::m_items;
uint64_t sdsl::mm::m_total_memory;
uint64_t* sdsl::mm::m_data;
std::ostream* sdsl::mm::m_out;
sdsl::util::stop_watch sdsl::mm::m_sw;
uint64_t sdsl::mm::m_granularity;
uint64_t sdsl::mm::m_pre_rtime;
uint64_t sdsl::mm::m_pre_max_mem;

sdsl::mm_initializer::mm_initializer()
{
    if (0 == nifty_counter++) {
        mm::m_total_memory = 0;
        mm::m_granularity = 0;
        mm::m_pre_rtime = 0;
        mm::m_pre_max_mem = 0;
        // initialize static members object here
        // mm::m_items.clear();
        mm::m_items = mm::tMVecItem();
        mm::m_data = NULL;
        mm::m_out = NULL;
    }
}
sdsl::mm_initializer::~mm_initializer()
{
    if (0 == --nifty_counter) {
        // clean up
    }
}

//! Namespace for the succinct data structure library
namespace sdsl
{


bool mm::map_hp()
{
#ifdef MAP_HUGETLB
    size_t hpgs= (m_total_memory+HUGE_LEN-1)/HUGE_LEN; // number of huge pages required to store the int_vectors
    m_data = (uint64_t*)mmap(NULL, hpgs*HUGE_LEN, HUGE_PROTECTION, HUGE_FLAGS, 0, 0);
    if (m_data == MAP_FAILED) {
        std::cout << "mmap was not successful" << std::endl;
        return false;
    }
    // map int_vectors
    uint64_t* addr = m_data;
    bool success = true;
    for (tMVecItem::const_iterator it=m_items.begin(); it!=m_items.end(); ++it) {
        success = success && it->second->map_hp(addr);
    }
    return success;
#else
    return false;
#endif
}

bool mm::unmap_hp()
{
#ifdef MAP_HUGETLB
    size_t hpgs= (m_total_memory+HUGE_LEN-1)/HUGE_LEN; // number of huge pages
    bool success = true;
    for (tMVecItem::const_iterator it=m_items.begin(); it!=m_items.end(); ++it) {
        success = success && it->second->unmap_hp();
    }
//		uint64_t* tmp_data = (uint64_t*)malloc(m_total_memory); // allocate memory for int_vectors
//		memcpy(tmp_data, m_data, len); // copy data from the mmapped region
    int ret = munmap((void*)m_data, hpgs*HUGE_LEN);
    if (ret == -1) {
        perror("Unmap failed");
        return false;
    }
    return success;
#else
    return true;
#endif
}

void mm::log_stream(std::ostream* out)
{
    m_out = out;
}

void mm::log_granularity(uint64_t granularity)
{
    m_granularity = granularity;
}

} // end namespace
