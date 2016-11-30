/*! \file ram_fs.hpp
 * \brief ram_fs.hpp
 * \author Simon Gog
 */
#ifndef INCLUDED_SDSL_RAM_FS
#define INCLUDED_SDSL_RAM_FS

#include "uintx_t.hpp"
#include "memory_tracking.hpp"
#include <string>
#include <map>
#include <vector>
#include <mutex>

namespace sdsl
{

class ram_fs;

//! ram_fs is a simple store for RAM-files.
/*!
 * (strings) to file content (content_type).
 */
class ram_fs
{
    public:
        typedef std::vector<char, track_allocator<char>> content_type;

    private:
        typedef std::map<std::string, content_type> mss_type;
        typedef std::map<int, std::string> mis_type;
        mss_type m_map;
        std::recursive_mutex m_rlock;
        mis_type m_fd_map;

        static ram_fs& the_ramfs() {
            static ram_fs fs;
            return fs;
        }
    public:
        //! Default construct
        ram_fs();
        static void store(const std::string& name, content_type data);
        //! Check if the file exists
        static bool exists(const std::string& name);
        //! Get the file size
        static size_t file_size(const std::string& name);

        //! Get the content
        static content_type& content(const std::string& name);
        //! Remove the file with key `name`
        static int remove(const std::string& name);
        //! Rename the file. Change key `old_filename` into `new_filename`.
        static int rename(const std::string old_filename, const std::string new_filename);

        //! Get fd for file
        static int open(const std::string& name);
        //! Get fd for file
        static int close(const int fd);
        //! Get the content with fd
        static content_type& content(const int fd);
        //! Get the content with fd
        static int truncate(const int fd,size_t new_size);
        //! Get the file size with fd_
        static size_t file_size(const int fd);
};

//! Determines if the given file is a RAM-file.
bool is_ram_file(const std::string& file);

//! Determines if the given file is a RAM-file.
bool is_ram_file(const int fd);

//! Returns the corresponding RAM-file name for file.
std::string ram_file_name(const std::string& file);

//! Returns for a RAM-file the corresponding disk file name
std::string disk_file_name(const std::string& file);

//! Remove a file.
int remove(const std::string& file);

//! Rename a file
int rename(const std::string& old_filename, const std::string& new_filename);

} // end namespace sdsl
#endif
