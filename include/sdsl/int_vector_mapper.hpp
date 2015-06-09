#ifndef SDSL_INT_VECTOR_MAPPER
#define SDSL_INT_VECTOR_MAPPER

#include "int_vector.hpp"
#include <sys/mman.h>
#include <fcntl.h>
#include <cstdio>
#include <ios>

namespace sdsl
{

template <uint8_t t_width = 0,std::ios_base::openmode t_mode = std::ios_base::out|std::ios_base::in>
class int_vector_mapper
{
        static_assert(t_width <= 64,
                      "int_vector_mapper: width must be at most 64 bits.");
    public:
        typedef typename int_vector<t_width>::difference_type difference_type;
        typedef typename int_vector<t_width>::value_type value_type;
        typedef typename int_vector<t_width>::size_type size_type;
        typedef typename int_vector<t_width>::int_width_type width_type;
    public:
        const size_type append_block_size = 1000000;
    private:
        uint8_t* m_mapped_data = nullptr;
        uint64_t m_file_size_bytes = 0;
        off_t m_data_offset = 0;
        int m_fd = -1;
        int_vector<t_width> m_wrapper;
        std::string m_file_name;
        bool m_delete_on_close;
    private:
        void mmap_file()
        {
            if (!(t_mode&std::ios_base::out)) {  // read only
                m_mapped_data = (uint8_t*)mmap(NULL,
                                               m_file_size_bytes,
                                               PROT_READ,
                                               MAP_SHARED,
                                               m_fd, 0);
            } else {
                m_mapped_data = (uint8_t*)mmap(NULL,
                                               m_file_size_bytes,
                                               PROT_READ|PROT_WRITE,
                                               MAP_SHARED,
                                               m_fd, 0);
            }
            if (m_mapped_data == MAP_FAILED) {
                std::string mmap_error
                    = std::string("int_vector_mapper: mmap error. ")
                      + std::string(strerror(errno));
                throw std::runtime_error(mmap_error);
            }
            auto ret = madvise(m_mapped_data, m_file_size_bytes, MADV_SEQUENTIAL);
            if (ret == -1) {
                perror("Error trying to hint sequential access");
            }
        }
    public:
        int_vector_mapper() = delete;
        int_vector_mapper(const int_vector_mapper&) = delete;
        int_vector_mapper& operator=(const int_vector_mapper&) = delete;
    public:
        ~int_vector_mapper()
        {
            if (m_mapped_data) {
                if (t_mode&std::ios_base::out) { // write was possible
                    if (m_data_offset) {
                        // update size in the on disk representation and
                        // truncate if necessary
                        uint64_t* size_in_file = (uint64_t*)m_mapped_data;
                        if (*size_in_file != m_wrapper.m_size) {
                            *size_in_file = m_wrapper.m_size;
                        }
                        if (t_width==0) {
                            // if size is variable and we map a sdsl vector
                            // we might have to update the stored width
                            uint8_t stored_width = m_mapped_data[8];
                            if (stored_width != m_wrapper.m_width) {
                                m_mapped_data[8] = m_wrapper.m_width;
                            }
                        }
                    }
                }

                munmap(m_mapped_data, m_file_size_bytes);

                if (t_mode&std::ios_base::out) {
                    // do we have to truncate?
                    size_type current_bit_size = m_wrapper.m_size;
                    size_type data_size_in_bytes = ((current_bit_size + 63) >> 6) << 3;
                    if (m_file_size_bytes != data_size_in_bytes + m_data_offset) {
                        int tret = ftruncate(m_fd, data_size_in_bytes + m_data_offset);
                        if (tret == -1) {
                            std::string truncate_error
                                = std::string("int_vector_mapper: truncate error. ")
                                  + std::string(strerror(errno));
                            throw std::runtime_error(truncate_error);
                        }
                    }
                }
            }
            if (m_fd != -1) {
                close(m_fd);
                if (m_delete_on_close) {
                    int ret_code = sdsl::remove(m_file_name);
                    if (ret_code != 0) {
                        std::cerr << "int_vector_mapper: error deleting file '"
                                  << m_file_name << "': " << ret_code << std::endl;
                    }
                }
            }
            m_wrapper.m_data = nullptr;
            m_wrapper.m_size = 0;
        }
        int_vector_mapper(int_vector_mapper&& ivm)
        {
            m_wrapper.m_data = ivm.m_wrapper.m_data;
            m_wrapper.m_size = ivm.m_wrapper.m_size;
            m_wrapper.width(ivm.m_wrapper.width());
            m_file_name = ivm.m_file_name;
            m_delete_on_close = ivm.m_delete_on_close;
            ivm.m_wrapper.m_data = nullptr;
            ivm.m_wrapper.m_size = 0;
            ivm.m_mapped_data = nullptr;
            ivm.m_fd = -1;
        }
        int_vector_mapper& operator=(int_vector_mapper&& ivm)
        {
            m_wrapper.m_data = ivm.m_wrapper.m_data;
            m_wrapper.m_size = ivm.m_wrapper.m_size;
            m_wrapper.width(ivm.m_wrapper.width());
            m_file_name = ivm.m_file_name;
            m_delete_on_close = ivm.m_delete_on_close;
            ivm.m_wrapper.m_data = nullptr;
            ivm.m_wrapper.m_size = 0;
            ivm.m_mapped_data = nullptr;
            ivm.m_fd = -1;
            return (*this);
        }
        int_vector_mapper(const std::string& key,const cache_config& config)
            : int_vector_mapper(cache_file_name(key, config)) {}
        int_vector_mapper(const std::string filename,
                          bool is_plain = false,
                          bool delete_on_close = false) :
            m_file_name(filename), m_delete_on_close(delete_on_close)
        {
            size_type size_in_bits = 0;
            uint8_t int_width = t_width;
            {
                std::ifstream f(filename);
                if (!f.is_open()) {
                    throw std::runtime_error(
                        "int_vector_mapper: file does not exist.");
                }
                if (!is_plain) {
                    int_vector<t_width>::read_header(size_in_bits, int_width, f);
                }
            }
            m_file_size_bytes = util::file_size(filename);

            if (!is_plain) {
                m_data_offset = t_width ? 8 : 9;
            } else {
                if (8 != t_width and 16 != t_width and 32 != t_width and 64
                    != t_width) {
                    throw std::runtime_error("int_vector_mapper: plain vector can "
                                             "only be of width 8, 16, 32, 64.");
                }
                size_in_bits = m_file_size_bytes * 8;
            }

            // open backend file depending on mode
            if (!(t_mode&std::ios_base::out))
                m_fd = open(filename.c_str(), O_RDONLY);
            else
                m_fd = open(filename.c_str(), O_RDWR);

            if (m_fd == -1) {
                std::string open_error
                    = std::string("int_vector_mapper: open file error.")
                      + std::string(strerror(errno));
                throw std::runtime_error(open_error);
            }

            // prepare for mmap
            size_type new_size_in_bytes = ((size_in_bits + 63) >> 6) << 3;
            m_file_size_bytes = new_size_in_bytes + m_data_offset;
            m_wrapper.width(int_width);
            mmap_file();
            m_wrapper.m_size = size_in_bits;
            m_wrapper.m_data = (uint64_t*)(m_mapped_data + m_data_offset);
        }
        std::string file_name() const { return m_file_name; }
        width_type width() const { return m_wrapper.width(); }
        void width(const uint8_t new_int_width)
        {
            static_assert(t_mode & std::ios_base::out,"int_vector_mapper: must be opened in in+out mode for 'width'");
            m_wrapper.width(new_int_width);
        }
        size_type size() const
        {
            return m_wrapper.size();
        }
        void bit_resize(const size_type bit_size)
        {
            static_assert(t_mode & std::ios_base::out,"int_vector_mapper: must be opened in in+out mode for 'bit_resize'");
            size_type new_size_in_bytes = ((bit_size + 63) >> 6) << 3;
            if (m_file_size_bytes != new_size_in_bytes + m_data_offset) {
                if (m_mapped_data) munmap(m_mapped_data, m_file_size_bytes);
                int tret = ftruncate(m_fd, new_size_in_bytes + m_data_offset);
                if (tret == -1) {
                    std::string truncate_error
                        = std::string("int_vector_mapper: truncate error. ")
                          + std::string(strerror(errno));
                    throw std::runtime_error(truncate_error);
                }
                m_file_size_bytes = new_size_in_bytes + m_data_offset;
            }

            // perform the actual mapping
            mmap_file();

            // update wrapper
            m_wrapper.m_data = (uint64_t*)(m_mapped_data + m_data_offset);
            m_wrapper.m_size = bit_size;
        }


        void resize(const size_type size)
        {
            static_assert(t_mode & std::ios_base::out,"int_vector_mapper: must be opened in in+out mode for 'resize'");
            size_type size_in_bits = size * width();
            bit_resize(size_in_bits);
        }

        auto begin() -> typename int_vector<t_width>::iterator {
            static_assert(t_mode & std::ios_base::out,"int_vector_mapper: must be opened in in+out mode for 'begin'");
            return m_wrapper.begin();
        }
        auto end() -> typename int_vector<t_width>::iterator {
            static_assert(t_mode & std::ios_base::out,"int_vector_mapper: must be opened in in+out mode for 'end'");
            return m_wrapper.end();
        }
        auto begin() const -> typename int_vector<t_width>::const_iterator {
            return m_wrapper.begin();
        }
        auto end() const -> typename int_vector<t_width>::const_iterator {
            return m_wrapper.end();
        }
        auto cbegin() const -> typename int_vector<t_width>::const_iterator {
            return m_wrapper.begin();
        }
        auto cend() const -> typename int_vector<t_width>::const_iterator {
            return m_wrapper.end();
        }
        auto operator[](const size_type& idx) const
        -> typename int_vector<t_width>::const_reference {
            return m_wrapper[idx];
        }
        auto operator[](const size_type& idx)
        -> typename int_vector<t_width>::reference {
            static_assert(t_mode & std::ios_base::out,"int_vector_mapper: must be opened in in+out mode for 'operator[]'");
            return m_wrapper[idx];
        }
        const uint64_t* data() const { return m_wrapper.data(); }
        uint64_t* data()
        {
            static_assert(t_mode & std::ios_base::out,"int_vector_mapper: must be opened in in+out mode for 'data'");
            return m_wrapper.data();
        }
        value_type get_int(size_type idx, const uint8_t len = 64) const
        {
            return m_wrapper.get_int(idx, len);
        }
        void set_int(size_type idx, value_type x, const uint8_t len = 64)
        {
            static_assert(t_mode & std::ios_base::out,"int_vector_mapper: must be opened in in+out mode for 'set_int'");
            m_wrapper.set_int(idx, x, len);
        }
        void push_back(value_type x)
        {
            static_assert(t_mode & std::ios_base::out,"int_vector_mapper: must be opened in in+out mode for 'push_back'");
            if (capacity() < size() + 1) {
                size_type old_size = m_wrapper.m_size;
                size_type size_in_bits = (size() + append_block_size) * width();
                bit_resize(size_in_bits);
                m_wrapper.m_size = old_size;
            }
            // update size in wrapper only
            m_wrapper.m_size += width();
            m_wrapper[size()-1] = x;
        }
        size_type capacity() const
        {
            size_t data_size_in_bits = 8 * (m_file_size_bytes - m_data_offset);
            return data_size_in_bits / width();
        }
        size_type bit_size() const
        {
            return m_wrapper.bit_size();
        }
        template<class container>
        bool operator==(const container& v) const
        {
            return std::equal(begin(), end(), v.begin());
        }
        bool operator==(const int_vector<t_width>& v) const
        {
            return m_wrapper == v;
        }
        bool operator==(const int_vector_mapper& v) const
        {
            return m_wrapper == v.m_wrapper;
        }
        template<class container>
        bool operator!=(const container& v) const
        {
            return !(*this==v);
        }
        void flip()
        {
            static_assert(t_mode & std::ios_base::out,"int_vector_mapper: must be opened in in+out mode for 'flip'");
            m_wrapper.flip();
        }
        bool empty() const
        {
            return m_wrapper.empty();
        }
};

template <uint8_t t_width = 0>
class temp_file_buffer
{
    private:
        static std::string tmp_file(const std::string& dir)
        {
            char tmp_file_name[1024] = {0};
            sprintf(tmp_file_name, "%s/tmp_mapper_file_XXXXXX.sdsl",dir.c_str());
            int fd = mkstemps(tmp_file_name,5);
            if (fd == -1) {
                throw std::runtime_error("could not create temporary file.");
            }
            close(fd);
            return std::string(tmp_file_name,strlen(tmp_file_name));
        }
    public:
        static int_vector_mapper<t_width> create()
        {
            auto file_name = tmp_file("/tmp");
            return create(file_name);
        }
        static int_vector_mapper<t_width> create(const cache_config& config)
        {
            auto file_name = tmp_file(config.dir);
            return create(file_name);
        }
        static int_vector_mapper<t_width> create(const std::string& file_name)
        {
            //write empty int_vector to init the file
            int_vector<t_width> tmp_vector;
            store_to_file(tmp_vector,file_name);
            return int_vector_mapper<t_width,std::ios_base::out|std::ios_base::in>(file_name,false,true);
        }
};

// creates emtpy int_vector<> that will not be deleted
template <uint8_t t_width = 0>
class write_out_buffer
{
    public:
        static int_vector_mapper<t_width> create(const std::string& key,cache_config& config)
        {
            auto file_name = cache_file_name(key,config);
            auto tmp = create(file_name);
            register_cache_file(key,config);
            return std::move(tmp);
        }
        static int_vector_mapper<t_width> create(const std::string& file_name)
        {
            //write empty int_vector to init the file
            int_vector<t_width> tmp_vector;
            store_to_file(tmp_vector,file_name);
            return int_vector_mapper<t_width,std::ios_base::out|std::ios_base::in>(file_name,false,false);
        }
};

template<std::ios_base::openmode t_mode = std::ios_base::out|std::ios_base::in>
using bit_vector_mapper = int_vector_mapper<1,t_mode>;

} // end of namespace

#endif
