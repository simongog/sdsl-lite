/*!\file memory_management.hpp
\brief memory_management.hpp contains two function for allocating and deallocating memory
\author Simon Gog
*/
#ifndef INCLUDED_SDSL_MEMORY_MANAGEMENT
#define INCLUDED_SDSL_MEMORY_MANAGEMENT

#include "uintx_t.hpp"
#include "config.hpp"
#include "bits.hpp"
#include "memory_tracking.hpp"
#include "ram_fs.hpp"


namespace sdsl
{

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


#ifndef MSVC_COMPILER
class hugepage_allocator
{
    private:
        uint8_t* m_base = nullptr;
        mm_block_t* m_first_block = nullptr;
        uint8_t* m_top = nullptr;
        size_t m_total_size = 0;
        std::multimap<size_t, mm_block_t*> m_free_large;
    private:
        size_t determine_available_hugepage_memory();
        void coalesce_block(mm_block_t* block);
        void split_block(mm_block_t* bptr, size_t size);
        uint8_t* hsbrk(size_t size);
        mm_block_t* new_block(size_t size);
        void remove_from_free_set(mm_block_t* block);
        void insert_into_free_set(mm_block_t* block);
        mm_block_t* find_free_block(size_t size_in_bytes);
        mm_block_t* last_block();
        void print_heap();
    public:
        void init(SDSL_UNUSED size_t size_in_bytes = 0)
        {
#ifdef MAP_HUGETLB
            if (size_in_bytes == 0) {
                size_in_bytes = determine_available_hugepage_memory();
            }

            m_total_size = size_in_bytes;
            m_base = (uint8_t*)mmap(nullptr, m_total_size,
                                    (PROT_READ | PROT_WRITE),
                                    (MAP_HUGETLB | MAP_ANONYMOUS | MAP_PRIVATE), 0, 0);
            if (m_base == MAP_FAILED) {
                throw std::system_error(ENOMEM, std::system_category(),
                                        "hugepage_allocator could not allocate hugepages");
            } else {
                // init the allocator
                m_top = m_base;
                m_first_block = (mm_block_t*)m_base;
            }
#else
            throw std::system_error(ENOMEM, std::system_category(),
                                    "hugepage_allocator: MAP_HUGETLB / hugepage support not available");
#endif
        }
        void* mm_realloc(void* ptr, size_t size);
        void* mm_alloc(size_t size_in_bytes);
        void mm_free(void* ptr);
        bool in_address_space(void* ptr)
        {
            // check if ptr is in the hugepage address space
            if (ptr == nullptr) {
                return true;
            }
            if (ptr >= m_base && ptr < m_top) {
                return true;
            }
            return false;
        }
        static hugepage_allocator& the_allocator()
        {
            static hugepage_allocator a;
            return a;
        }
};
#endif

class memory_manager
{
    private:
        bool hugepages = false;
    private:
        static memory_manager& the_manager()
        {
            static memory_manager m;
            return m;
        }
    public:
        static uint64_t* alloc_mem(size_t size_in_bytes)
        {
#ifndef MSVC_COMPILER
            auto& m = the_manager();
            if (m.hugepages) {
                return (uint64_t*)hugepage_allocator::the_allocator().mm_alloc(size_in_bytes);
            }
#endif
            return (uint64_t*)calloc(size_in_bytes, 1);
        }
        static void free_mem(uint64_t* ptr)
        {
#ifndef MSVC_COMPILER
            auto& m = the_manager();
            if (m.hugepages and hugepage_allocator::the_allocator().in_address_space(ptr)) {
                hugepage_allocator::the_allocator().mm_free(ptr);
                return;
            }
#endif
            std::free(ptr);
        }
        static uint64_t* realloc_mem(uint64_t* ptr, size_t size)
        {
#ifndef MSVC_COMPILER
            auto& m = the_manager();
            if (m.hugepages and hugepage_allocator::the_allocator().in_address_space(ptr)) {
                return (uint64_t*)hugepage_allocator::the_allocator().mm_realloc(ptr, size);
            }
#endif
            uint64_t* temp = (uint64_t*)realloc(ptr, size);
            if (temp == NULL) {
                throw std::bad_alloc();
            }
            return temp;
        }
    public:
        static void use_hugepages(size_t bytes = 0)
        {
#ifndef MSVC_COMPILER
            auto& m = the_manager();
            hugepage_allocator::the_allocator().init(bytes);
            m.hugepages = true;
#else
            throw std::runtime_error("hugepages not support on MSVC_COMPILER");
#endif
        }
        template<class t_vec>
        static void resize(t_vec& v, const typename t_vec::size_type size)
        {
            uint64_t old_size_in_bytes = ((v.m_size + 63) >> 6) << 3;
            uint64_t new_size_in_bytes = ((size + 63) >> 6) << 3;
            bool do_realloc = old_size_in_bytes != new_size_in_bytes;
            v.m_size = size;
            if (do_realloc || v.m_data == nullptr) {
                // Note that we allocate 8 additional bytes if m_size % 64 == 0.
                // We need this padding since rank data structures do a memory
                // access to this padding to answer rank(size()) if size()%64 ==0.
                // Note that this padding is not counted in the serialize method!
                size_t allocated_bytes = (size_t)(((size + 64) >> 6) << 3);
                v.m_data = memory_manager::realloc_mem(v.m_data, allocated_bytes);
                if (allocated_bytes != 0 && v.m_data == nullptr) {
                    throw std::bad_alloc();
                }
                // update and fill with 0s
                if (v.bit_size() < v.capacity()) {
                    uint8_t len = (uint8_t)(v.capacity() - v.bit_size());
                    uint8_t in_word_offset = (uint8_t)(v.bit_size() & 0x3F);
                    bits::write_int(v.m_data + (v.bit_size() >> 6), 0, in_word_offset, len);
                }
                if (((v.m_size) % 64) == 0) {  // initialize unreachable bits with 0
                    v.m_data[v.m_size / 64] = 0;
                }

                // update stats
                if (do_realloc) {
                    memory_monitor::record((int64_t)new_size_in_bytes - (int64_t)old_size_in_bytes);
                }
            }
        }
        template<class t_vec>
        static void clear(t_vec& v)
        {
            int64_t size_in_bytes = ((v.m_size + 63) >> 6) << 3;
            // remove mem
            memory_manager::free_mem(v.m_data);
            v.m_data = nullptr;

            // update stats
            if (size_in_bytes) {
                memory_monitor::record(size_in_bytes*-1);
            }
        }

        static int open_file_for_mmap(std::string& filename, std::ios_base::openmode mode) {
            if( is_ram_file(filename) ) {
                return ram_fs::open(filename);
            }
#ifdef MSVC_COMPILER
            int fd = -1;
            if (!(mode&std::ios_base::out)) _sopen_s(&fd,filename.c_str(), _O_BINARY| _O_RDONLY, _SH_DENYNO, _S_IREAD);
            else _sopen_s(&fd, filename.c_str(), _O_BINARY | _O_RDWR, _SH_DENYNO, _S_IREAD | _S_IWRITE);
            return fd;
#else
            if (!(mode&std::ios_base::out)) return open(filename.c_str(), O_RDONLY);
            else return open(filename.c_str(), O_RDWR);
#endif
            return -1;
        }

        static void* mmap_file(int fd,uint64_t file_size, std::ios_base::openmode mode) {
            if (file_size==0){
                std::cout<<"file_size=0"<<std::endl;
                return nullptr;
            }
            if( is_ram_file(fd) ) {
                if( ram_fs::file_size(fd) < file_size) return nullptr;
                auto& file_content = ram_fs::content(fd);
                return file_content.data();
            }
            memory_monitor::record(file_size);
#ifdef MSVC_COMPILER
            HANDLE fh = (HANDLE)_get_osfhandle(fd);
            if (fh == INVALID_HANDLE_VALUE) {
                return nullptr;
            }
            HANDLE fm;
            if (!(mode&std::ios_base::out)) { // read only?
                fm = CreateFileMapping(fh, NULL, PAGE_READONLY, 0, 0, NULL);
            } else fm = CreateFileMapping(fh, NULL, PAGE_READWRITE, 0, 0, NULL);
            if (fm == NULL) {
                return nullptr;
            }
            void* map = nullptr;
            if (!(mode&std::ios_base::out)) { // read only?
                map = MapViewOfFile(fm, FILE_MAP_READ, 0, 0, file_size);
            } else map = MapViewOfFile(fm, FILE_MAP_WRITE | FILE_MAP_READ, 0, 0, file_size);
            // we can close the file handle before we unmap the view: (see UnmapViewOfFile Doc)
            // Although an application may close the file handle used to create a file mapping object, 
            // the system holds the corresponding file open until the last view of the file is unmapped. 
            // Files for which the last view has not yet been unmapped are held open with no sharing restrictions.
            CloseHandle(fm);
            return map;
#else
            void* map = nullptr;
            if (!(mode&std::ios_base::out)) map = mmap(NULL,file_size,PROT_READ,MAP_SHARED,fd, 0);
            else map = mmap(NULL,file_size,PROT_READ | PROT_WRITE,MAP_SHARED,fd, 0);
            if(map == MAP_FAILED) map = nullptr; // unify windows and unix error behaviour
            return map;
#endif
            return nullptr;
        }

        static int mem_unmap(int fd,void* addr,const uint64_t size) {
            if ( addr == nullptr ) {
                return 0;
            }
            if( is_ram_file(fd) ) {
                return 0;
            }
            memory_monitor::record(-((int64_t)size));
#ifdef MSVC_COMPILER
            if (UnmapViewOfFile(addr)) return 0;
            return -1;
#else
            return munmap(addr, size);
#endif
            return -1;
        }

        static int close_file_for_mmap(int fd) {
            if( is_ram_file(fd) ) {
                return ram_fs::close(fd);
            }
#ifdef MSVC_COMPILER
            return _close(fd);
#else
            return close(fd);
#endif
            return -1;
        }

        static int truncate_file_mmap(int fd,const uint64_t new_size) {
            if( is_ram_file(fd) ) {
                return ram_fs::truncate(fd,new_size);
            }
#ifdef MSVC_COMPILER
            auto ret = _chsize_s(fd,new_size);
            if(ret != 0) ret = -1;
            return ret;
#else
            return ftruncate(fd,new_size);
#endif
            return -1;
        }

};

} // end namespace

#endif
