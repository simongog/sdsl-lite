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
#include <stack>
#include "config.hpp"

namespace sdsl
{

class memory_monitor;

template<format_type F>
void write_mem_log(std::ostream& out,const memory_monitor& m);

class memory_monitor
{
    public:
        using timer = std::chrono::high_resolution_clock;
        struct mm_alloc {
            timer::time_point timestamp;
            int64_t usage;
            mm_alloc(timer::time_point t,int64_t u) : timestamp(t) , usage(u) {};
        };
        struct mm_event {
            std::string name;
            std::vector<mm_alloc> allocations;
            mm_event(std::string n,int64_t usage) : name(n) {
                allocations.emplace_back(timer::now(),usage);
            };
            bool operator< (const mm_event& a) const {
                if (a.allocations.size() && this->allocations.size()) {
                    if (this->allocations[0].timestamp == a.allocations[0].timestamp) {
                        return this->allocations.back().timestamp < a.allocations.back().timestamp;
                    } else {
                        return this->allocations[0].timestamp < a.allocations[0].timestamp;
                    }
                }
                return true;
            }
        };
        struct mm_event_proxy {
            bool add;
            timer::time_point created;
            mm_event_proxy(const std::string& name,int64_t usage,bool a) : add(a) {
                if (add) {
                    auto& m = the_monitor();
                    std::lock_guard<util::spin_lock> lock(m.spinlock);
                    m.event_stack.emplace(name,usage);
                }
            }
            ~mm_event_proxy() {
                if (add) {
                    auto& m = the_monitor();
                    std::lock_guard<util::spin_lock> lock(m.spinlock);
                    auto& cur = m.event_stack.top();
                    auto cur_time = timer::now();
                    cur.allocations.emplace_back(cur_time,m.current_usage);
                    m.completed_events.emplace_back(std::move(cur));
                    m.event_stack.pop();
                    // add a point to the new "top" with the same memory
                    // as before but just ahead in time
                    if (! m.event_stack.empty()) {
                        if (m.event_stack.top().allocations.size()) {
                            auto last_usage = m.event_stack.top().allocations.back().usage;
                            m.event_stack.top().allocations.emplace_back(cur_time,last_usage);
                        }
                    }
                }
            }
        };
        std::chrono::milliseconds log_granularity = std::chrono::milliseconds(20);
        int64_t current_usage = 0;
        bool track_usage = false;
        std::vector<mm_event> completed_events;
        std::stack<mm_event> event_stack;
        timer::time_point start_log;
        timer::time_point last_event;
        util::spin_lock spinlock;
    private:
        // disable construction of the object
        memory_monitor() {};
        ~memory_monitor() {
            if (track_usage) {
                stop();
            }
        }
        memory_monitor(const memory_monitor&) = delete;
        memory_monitor& operator=(const memory_monitor&) = delete;
    private:
        static memory_monitor& the_monitor() {
            static memory_monitor m;
            return m;
        }
    public:
        static void granularity(std::chrono::milliseconds ms) {
            auto& m = the_monitor();
            m.log_granularity = ms;
        }
        static int64_t peak() {
            auto& m = the_monitor();
            int64_t max = 0;
            for (auto events : m.completed_events) {
                for (auto alloc : events.allocations) {
                    if (max < alloc.usage) {
                        max = alloc.usage;
                    }
                }
            }
            return max;
        }

        static void start() {
            auto& m = the_monitor();
            m.track_usage = true;
            // clear if there is something there
            if (m.completed_events.size()) {
                m.completed_events.clear();
            }
            while (m.event_stack.size()) {
                m.event_stack.pop();
            }
            m.start_log = timer::now();
            m.current_usage = 0;
            m.last_event = m.start_log;
            m.event_stack.emplace("unknown",0);
        }
        static void stop() {
            auto& m = the_monitor();
            while (! m.event_stack.empty()) {
                m.completed_events.emplace_back(std::move(m.event_stack.top()));
                m.event_stack.pop();
            }
            m.track_usage = false;
        }
        static void record(int64_t delta) {
            auto& m = the_monitor();
            if (m.track_usage) {
                std::lock_guard<util::spin_lock> lock(m.spinlock);
                auto cur = timer::now();
                if (m.last_event + m.log_granularity < cur) {
                    m.event_stack.top().allocations.emplace_back(cur,m.current_usage);
                    m.current_usage = m.current_usage + delta;
                    m.event_stack.top().allocations.emplace_back(cur,m.current_usage);
                    m.last_event = cur;
                } else {
                    if (m.event_stack.top().allocations.size()) {
                        m.current_usage = m.current_usage + delta;
                        m.event_stack.top().allocations.back().usage = m.current_usage;
                        m.event_stack.top().allocations.back().timestamp = cur;
                    }
                }
            }
        }
        static mm_event_proxy event(const std::string& name) {
            auto& m = the_monitor();
            if (m.track_usage) {
                return mm_event_proxy(name,m.current_usage,true);
            }
            return mm_event_proxy(name,m.current_usage,false);
        }
        template<format_type F>
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
        size_t m_total_size = 0;
        std::multimap<size_t,mm_block_t*> m_free_large;
    private:
        size_t determine_available_hugepage_memory();
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
        void init(SDSL_UNUSED size_t size_in_bytes = 0) {
#ifdef MAP_HUGETLB
            if (size_in_bytes == 0) {
                size_in_bytes = determine_available_hugepage_memory();
            }

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
        bool in_address_space(void* ptr) {
            // check if ptr is in the hugepage address space
            if (ptr == nullptr) {
                return true;
            }
            if (ptr >= m_base && ptr < m_top) {
                return true;
            }
            return false;
        }
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
            if (m.hugepages and hugepage_allocator::the_allocator().in_address_space(ptr)) {
                hugepage_allocator::the_allocator().mm_free(ptr);
            } else {
                std::free(ptr);
            }
        }
        static uint64_t* realloc_mem(uint64_t* ptr,size_t size) {
            auto& m = the_manager();
            if (m.hugepages and hugepage_allocator::the_allocator().in_address_space(ptr)) {
                return (uint64_t*) hugepage_allocator::the_allocator().mm_realloc(ptr,size);
            } else {
                return (uint64_t*) realloc(ptr,size);
            }
        }
    public:
        static void use_hugepages(size_t bytes = 0) {
            auto& m = the_manager();
            hugepage_allocator::the_allocator().init(bytes);
            m.hugepages = true;
        }
        template<class t_vec>
        static void resize(t_vec& v, const typename t_vec::size_type size) {
            uint64_t old_size_in_bytes = ((v.m_size+63)>>6)<<3;
            uint64_t new_size_in_bytes = ((size+63)>>6)<<3;
            bool do_realloc = old_size_in_bytes != new_size_in_bytes;
            v.m_size = size;
            if (do_realloc || v.m_data == nullptr) {
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
                    memory_monitor::record((int64_t)new_size_in_bytes-(int64_t)old_size_in_bytes);
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
