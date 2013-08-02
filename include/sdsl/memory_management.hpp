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

namespace sdsl
{

class mm_item_base
{
    public:
        mm_item_base() {};
        virtual bool map_hp(uint64_t*&) {
            return false;
        };// const = 0;
        virtual bool unmap_hp() {
            return false;
        };// const = 0;
        virtual ~mm_item_base() {};
        virtual uint64_t size() {
            return 0;
        };
};

template<class int_vec_t>
class mm_item : public mm_item_base
{
    private:
        int_vec_t* m_v;
    public:
        explicit mm_item(int_vec_t* v):m_v(v) {}
        ~mm_item() { }

        //! Map content of int_vector to a hugepage starting at address addr
        /*!
         *  Details: The content of the corresponding int_vector of mm_item
         *           is copied into a hugepage starting at address addr.
         *           The string position in the hugepage is then increased
         *           by the number of bytes used by the int_vector.
         *           So addr is the new starting address for the next
         *           mm_item which have to be mapped.
         */
        bool map_hp(uint64_t*& addr) {
            uint64_t len = size();
            if (m_v->m_data != nullptr) {
                memcpy((char*)addr, m_v->m_data, len); // copy old data
                free(m_v->m_data);
                m_v->m_data = addr;
                addr += (len/8);
            }
            return true;
        }

        //!
        bool unmap_hp() {
            uint64_t len = size();
            if (util::verbose) {
                std::cerr<<"unmap int_vector of size "<< len <<std::endl;
                std::cerr<<"m_data="<<m_v->m_data<<std::endl;
            }
            uint64_t* tmp_data = (uint64_t*)malloc(len); // allocate memory for m_data
            memcpy(tmp_data, m_v->m_data, len); // copy data from the mmapped region
            m_v->m_data = tmp_data;
            return true;
        }

        uint64_t size() {
            return (((m_v->bit_size()+63)>>6)<<3);
        }
};

// initialization helper for class mm
class mm_initializer
{
    public:
        mm_initializer();
        ~mm_initializer();
};

} // end namespace sdsl

static sdsl::mm_initializer init_mm;

namespace sdsl
{

// memory management class
class mm
{
        friend class mm_initializer;
        typedef std::map<uint64_t, mm_item_base*> tMVecItem;
        static tMVecItem m_items;
        static uint64_t m_total_memory;
        static uint64_t* m_data;
        static std::ostream* m_out;
        static util::stop_watch m_sw;
        static uint64_t m_granularity;
        static uint64_t m_pre_max_mem;
        static uint64_t m_pre_rtime;
        static util::spin_lock m_spinlock;

    public:
        mm();

        template<class int_vec_t>
        static void add(int_vec_t* v, bool moved=false) {
            std::lock_guard<util::spin_lock> lock(m_spinlock);
            if (mm::m_items.find((uint64_t)v) == mm::m_items.end()) {
                mm_item_base* item = new mm_item<int_vec_t>(v);
                if (false and util::verbose) {
                    std::cout << "mm::add: add vector " << v << std::endl;
                    std::cout.flush();
                }
                mm::m_items[(uint64_t)v] = item;
                if (!moved and item->size()) {
                    log("");
                    m_total_memory += item->size(); // add space
                    log("");
                }
            } else {
                if (false and util::verbose) std::cout << "mm::add: mm_item is already in the set" << std::endl;
            }
        }

        template<class int_vec_t>
        static void realloc(int_vec_t& v, const typename int_vec_t::size_type size) {
            bool do_realloc = ((size+63)>>6) != ((v.m_size+63)>>6);
            uint64_t old_size = ((v.m_size+63)>>6)<<3;
            v.m_size = size;                         // set new size
            // special case: bitvector of size 0
            if (do_realloc or v.m_data==nullptr) { // or (t_width==1 and m_size==0) ) {
                uint64_t* data = nullptr;
                // Note that we allocate 8 additional bytes if m_size % 64 == 0.
                // We need this padding since rank data structures do a memory
                // access to this padding to answer rank(size()) if size()%64 ==0.
                // Note that this padding is not counted in the serialize method!
                data = (uint64_t*)::realloc(v.m_data, (((v.m_size+64)>>6)<<3)); // if m_data == nullptr realloc
                // Method realloc is equivalent to malloc if m_data == nullptr.
                // If size is zero and ptr is not nullptr, a new, minimum sized object is allocated and the original object is freed.
                // The allocated memory is aligned such that it can be used for any data type, including AltiVec- and SSE-related types.
                v.m_data = data;
                // initialize unreachable bits to 0
                if (v.bit_size() < v.capacity()) {   //m_size>0
                    bits::write_int(v.m_data+(v.bit_size()>>6), 0, v.bit_size()&0x3F, v.capacity() - v.bit_size());
                }
                if (((v.m_size) % 64) == 0) {  // initialize unreachable bits with 0
                    v.m_data[v.m_size/64] = 0;
                }
            }
            if (old_size != ((v.m_size+63)>>6)<<3) {
                log("");
                {
                    std::lock_guard<util::spin_lock> lock(m_spinlock);
                    m_total_memory -= old_size; // subtract old space
                    m_total_memory += ((v.m_size+63)>>6)<<3; // add new space
                }
                log("");
            }
        }

        template<class int_vec_t>
        static void remove(int_vec_t* v) {
            std::lock_guard<util::spin_lock> lock(m_spinlock);
            if (mm::m_items.find((uint64_t)v) != mm::m_items.end()) {
                if (false and util::verbose) {
                    std::cout << "mm:remove: remove vector " << v << std::endl;
                };
                mm_item_base* item = m_items[(uint64_t)v];
                if (item->size()) {
                    log("");
                    m_total_memory -= item->size(); // delete space
                    log("");
                }
                mm::m_items.erase((uint64_t)v);
                delete item;
            } else {
                if (false and util::verbose) {
                    std::cout << "mm:remove: mm_item is not in the set" << std::endl;
                };
            }
        }

        static void log_stream(std::ostream* out);

        static void log_granularity(uint64_t granularity);

        static void log(const std::string& msg) {
            if (m_out != nullptr) {
                m_sw.stop();
                uint64_t log_time = m_sw.abs_real_time();
                if (log_time >= m_pre_rtime + m_granularity
                    or msg.size() > 0) {
                    (*m_out) << m_pre_rtime << ";" << m_sw.abs_user_time()  << ";"
                             << m_sw.abs_sys_time()  << ";" << m_pre_max_mem << ";"
                             << "" << std::endl;
                    if (msg.size() > 0) {  // output if msg is set

                        (*m_out) << log_time << ";" << m_sw.abs_user_time()  << ";"
                                 << m_sw.abs_sys_time()  << ";" << m_total_memory << ";"
                                 << msg << std::endl;
                    }

                    m_pre_max_mem = m_total_memory; // reset memory
                    m_pre_rtime = log_time;
                } else {
                    m_pre_max_mem = std::max(m_pre_max_mem, m_total_memory);
                }
            }
        }

        static bool map_hp();
        static bool unmap_hp();
};

} // end namespace

#endif
