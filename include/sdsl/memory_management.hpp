/*!\file memory_management.hpp
   \brief memory_management.hpp contains two function for allocating and deallocating memory
   \author Simon Gog
*/
#ifndef INCLUDED_SDSL_MEMORY_MANAGEMENT
#define INCLUDED_SDSL_MEMORY_MANAGEMENT

#include "uintx_t.hpp"
#include "util.hpp"
#include <map>
#include <iostream>
#include <cstdlib>
#include <mutex>
#include <chrono>
#include <cstring>
#include <set>
#include <cstddef>


namespace sdsl
{

class memory_monitor;

enum memformat_type {JSON, CSV, HTML};

template<memformat_type F>
void write_mem_log(std::ostream& out,const memory_monitor& m);

class memory_monitor
{
    public:
        using timer = std::chrono::high_resolution_clock;
        struct mm_event {
            timer::time_point timestamp;
            int64_t usage;
            mm_event(timer::time_point t,int64_t u) : timestamp(t) , usage(u) {};
        };
        std::chrono::milliseconds log_granularity = std::chrono::milliseconds(20);
        int64_t current_usage = 0;
        int64_t peak_usage = 0;
        bool track_usage = false;
        std::vector<mm_event> mem_events;
        std::vector<std::tuple<timer::time_point,int64_t,std::string>> events;
        util::spin_lock spinlock;
    private:
        // disable construction of the object
        memory_monitor() {};
        memory_monitor(const memory_monitor&) = delete;
        memory_monitor& operator=(const memory_monitor&) = delete;
    private:
        static memory_monitor& the_monitor() {
            static memory_monitor m;
            return m;
        }
        static mm_event& last_event() {
            auto& m = the_monitor();
            if (!m.mem_events.size()) {
                m.mem_events.emplace_back(timer::now(),(int64_t)0);
            }
            return m.mem_events.back(); // empty event
        }
    public:
        static void granularity(std::chrono::milliseconds ms) {
            auto& m = the_monitor();
            m.log_granularity = ms;
        }
        static void start() {
            auto& m = the_monitor();
            m.track_usage = true;
            event("start mem_monitor");
        }
        static void stop() {
            auto& m = the_monitor();
            event("stop mem_monitor");
            m.mem_events.emplace_back(timer::now(),m.current_usage); // final event
            m.track_usage = false;
        }
        static void record(int64_t delta) {
            auto& m = the_monitor();
            if (m.track_usage) {
                std::lock_guard<util::spin_lock> lock(m.spinlock);
                m.current_usage = (int64_t)((int64_t)m.current_usage + delta);
                m.peak_usage = std::max(m.current_usage,m.peak_usage);
                auto cur = timer::now();
                if (last_event().timestamp + m.log_granularity < cur) {

                    m.mem_events.emplace_back(cur,m.current_usage);
                }
            }
        }
        static void event(const std::string& name) {
            auto& m = the_monitor();
            if (m.track_usage) {
                std::lock_guard<util::spin_lock> lock(m.spinlock);
                m.events.emplace_back(timer::now(),m.current_usage,name);
            }
        }
        template<memformat_type F>
        static void write_memory_log(std::ostream& out) {
            write_mem_log<F>(out,the_monitor());
        }
};

#pragma pack(push, 1)
typedef struct mm_block {
    size_t size;
    struct mm_block* next;
    struct mm_block* prev;
} mm_block_t;

typedef struct bfoot {
    size_t size;
} mm_block_foot_t;
#pragma pack(pop)

#include <sys/mman.h>

class hugepage_allocator
{
    private:
        uint8_t* m_base = nullptr;
        mm_block_t* m_first_block = nullptr;
        uint8_t* m_top = nullptr;
        size_t m_small_threshold = 4096;
        size_t m_total_size = 0;
        std::multimap<size_t,mm_block_t*> m_free_large;
    private:
        void coalesce_block(mm_block_t* block);
        void split_block(mm_block_t* bptr,size_t size);
        uint8_t* hsbrk(size_t size);
        mm_block_t* new_block(size_t size);
        void remove_from_free_set(mm_block_t* block);
        void insert_into_free_set(mm_block_t* block);
        mm_block_t* find_free_block(size_t size_in_bytes);
        mm_block_t* last_block();
        void print_heap();
    public:
        void init(SDSL_UNUSED size_t size_in_bytes) {
#ifdef MAP_HUGETLB
            m_total_size = size_in_bytes;
            m_base = (uint8_t*) mmap(nullptr, m_total_size,
                                     (PROT_READ | PROT_WRITE),
                                     (MAP_HUGETLB | MAP_ANONYMOUS | MAP_PRIVATE), 0, 0);
            if (m_base == MAP_FAILED) {
                throw std::system_error(ENOMEM,std::system_category(),
                                        "hugepage_allocator could not allocate hugepages");
            } else {
                // init the allocator
                m_top = m_base;
                m_first_block = (mm_block_t*) m_base;
            }
#else
            throw std::system_error(ENOMEM,std::system_category(),
                                    "hugepage_allocator: MAP_HUGETLB / hugepage support not available");
#endif
        }
        void* mm_realloc(void* ptr, size_t size);
        void* mm_alloc(size_t size_in_bytes);
        void mm_free(void* ptr);
        static hugepage_allocator& the_allocator() {
            static hugepage_allocator a;
            return a;
        }
};

class memory_manager
{
    private:
        bool hugepages = false;
    private:
        static memory_manager& the_manager() {
            static memory_manager m;
            return m;
        }
    public:
        static uint64_t* alloc_mem(size_t size_in_bytes) {
            auto& m = the_manager();
            if (m.hugepages) {
                return (uint64_t*) hugepage_allocator::the_allocator().mm_alloc(size_in_bytes);
            } else {
                return (uint64_t*) calloc(size_in_bytes,1);
            }
        }
        static void free_mem(uint64_t* ptr) {
            auto& m = the_manager();
            if (m.hugepages) {
                hugepage_allocator::the_allocator().mm_free(ptr);
            } else {
                std::free(ptr);
            }
        }
        static uint64_t* realloc_mem(uint64_t* ptr,size_t size) {
            auto& m = the_manager();
            if (m.hugepages) {
                return (uint64_t*) hugepage_allocator::the_allocator().mm_realloc(ptr,size);
            } else {
                return (uint64_t*) realloc(ptr,size);
            }
        }
    public:
        static void use_hugepages(size_t bytes) {
            auto& m = the_manager();
            m.hugepages = true;
            hugepage_allocator::the_allocator().init(bytes);
        }
        template<class t_vec>
        static void resize(t_vec& v, const typename t_vec::size_type size) {
            int64_t old_size_in_bytes = ((v.m_size+63)>>6)<<3;
            int64_t new_size_in_bytes = ((size+63)>>6)<<3;
            bool do_realloc = old_size_in_bytes != new_size_in_bytes;
            v.m_size = size;
            if (do_realloc || new_size_in_bytes == 0) {
                // Note that we allocate 8 additional bytes if m_size % 64 == 0.
                // We need this padding since rank data structures do a memory
                // access to this padding to answer rank(size()) if size()%64 ==0.
                // Note that this padding is not counted in the serialize method!
                size_t allocated_bytes = (((size+64)>>6)<<3);
                v.m_data = memory_manager::realloc_mem(v.m_data,allocated_bytes);
                if (allocated_bytes != 0 && v.m_data == nullptr) {
                    throw std::bad_alloc();
                }
                // update and fill with 0s
                if (v.bit_size() < v.capacity()) {
                    bits::write_int(v.m_data+(v.bit_size()>>6), 0, v.bit_size()&0x3F, v.capacity() - v.bit_size());
                }
                if (((v.m_size) % 64) == 0) {  // initialize unreachable bits with 0
                    v.m_data[v.m_size/64] = 0;
                }

                // update stats
                if (do_realloc) {
                    memory_monitor::record(new_size_in_bytes-old_size_in_bytes);
                }
            }
        }
        template<class t_vec>
        static void clear(t_vec& v) {
            int64_t size_in_bytes = ((v.m_size+63)>>6)<<3;
            // remove mem
            memory_manager::free_mem(v.m_data);
            v.m_data = nullptr;

            // update stats
            if (size_in_bytes) {
                memory_monitor::record(size_in_bytes*-1);
            }
        }
};


} // end namespace

#endif
