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
            uint64_t usage;
            mm_event(timer::time_point t,uint64_t u) : timestamp(t) , usage(u) {};
        };
        std::chrono::milliseconds log_granularity = std::chrono::milliseconds(50);
        uint64_t current_usage = 0;
        uint64_t peak_usage = 0;
        bool track_usage = false;
        std::vector<mm_event> mem_events;
        std::vector<std::pair<timer::time_point,std::string>> events;
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
                m.mem_events.emplace_back(timer::now(),(uint64_t)0);
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
                m.current_usage = (uint64_t)((int64_t)m.current_usage + delta);
                m.peak_usage = std::max(m.current_usage,m.peak_usage);
                auto cur = timer::now();
                if (last_event().timestamp + m.log_granularity > cur) {
                    m.mem_events.emplace_back(cur,m.current_usage);
                }
            }
        }
        static void event(const std::string& name) {
            auto& m = the_monitor();
            if (m.track_usage) {
                std::lock_guard<util::spin_lock> lock(m.spinlock);
                m.events.emplace_back(timer::now(),name);
            }
        }
        template<memformat_type F>
        static void write_memory_log(std::ostream& out) {
            write_mem_log<F>(out,the_monitor());
        }
};


#include <sys/mman.h>

class hugepage_allocator
{
    private:
        uint64_t* m_memory = nullptr;
        size_t m_mem_size = 0;
    private:
        static hugepage_allocator& the_allocator() {
            static hugepage_allocator a;
            return a;
        }
    public:
        static void init(size_t size_in_bytes) {
            auto& a = the_allocator();
            a.m_mem_size = size_in_bytes;
            a.m_memory = (uint64_t*) mmap(nullptr, size_in_bytes,
                                          (PROT_READ | PROT_WRITE),
                                          (MAP_HUGETLB | MAP_ANONYMOUS | MAP_PRIVATE), 0, 0);
            if (a.m_memory == MAP_FAILED) {
                throw std::bad_alloc();
            }
        }
        static uint64_t* alloc(size_t size_in_bytes) {
            return (uint64_t*) malloc(size_in_bytes);
        }
        static void free(uint64_t* ptr) {
            free(ptr);
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
        static uint64_t* alloc_mem(size_t size_in_bytes) {
            auto& m = the_manager();
            if (m.hugepages) {
                return hugepage_allocator::alloc(size_in_bytes);
            } else {
                return (uint64_t*) calloc(size_in_bytes,1);
            }
        }
        static void free_mem(uint64_t* ptr) {
            auto& m = the_manager();
            if (m.hugepages) {
                hugepage_allocator::free(ptr);
            } else {
                std::free(ptr);
            }
        }
    public:
        static void use_hugepages(size_t bytes) {
            auto& m = the_manager();
            m.hugepages = true;
            hugepage_allocator::init(bytes);
        }
        template<class t_vec>
        static void resize(t_vec& v, const typename t_vec::size_type size) {
            int64_t old_size_in_bytes = ((v.m_size+63)>>6)<<3;
            int64_t new_size_in_bytes = ((size+63)>>6)<<3;
            //std::cout << "resize(" << old_size_in_bytes << " , " << new_size_in_bytes << ")\n";
            bool do_realloc = old_size_in_bytes != new_size_in_bytes;
            if (do_realloc || new_size_in_bytes == 0) {
                // Note that we allocate 8 additional bytes if m_size % 64 == 0.
                // We need this padding since rank data structures do a memory
                // access to this padding to answer rank(size()) if size()%64 ==0.
                // Note that this padding is not counted in the serialize method!
                size_t allocated_bytes = (((size+64)>>6)<<3);
                uint64_t* data = memory_manager::alloc_mem(allocated_bytes);
                if (allocated_bytes != 0 && data == nullptr) {
                    throw std::bad_alloc();
                }
                // copy and update
                std::memcpy(data, v.m_data, std::min(old_size_in_bytes,new_size_in_bytes));
                memory_manager::free_mem(v.m_data);
                v.m_data = data;
                v.m_size = size;

                // update stats
                memory_monitor::record(new_size_in_bytes-old_size_in_bytes);
            }
            //std::cout << "done(" << old_size_in_bytes << " , " << new_size_in_bytes << ")\n";
        }
        template<class t_vec>
        static void clear(t_vec& v) {
            int64_t size_in_bytes = ((v.m_size+63)>>6)<<3;

            // remove mem
            memory_manager::free_mem(v.m_data);
            v.m_data = nullptr;

            // update stats
            memory_monitor::record(size_in_bytes*-1);
        }
};



} // end namespace

#endif
