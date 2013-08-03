/*! \file ram_fs.hpp
 * \brief ram_fs.hpp
 * \author Simon Gog
 */
#ifndef INCLUDED_SDSL_RAM_FS
#define INCLUDED_SDSL_RAM_FS

#include "uintx_t.hpp"
#include <string>
#include <map>
#include <vector>
#include <mutex>

namespace sdsl
{

class ram_fs_initializer
{
    public:
        ram_fs_initializer();
        ~ram_fs_initializer();
};

} // end namespace sdsl


static sdsl::ram_fs_initializer init_ram_fs;

namespace sdsl
{


//! ram_fs is a simple store for RAM-files.
/*!
 * Simple key-value store which maps file names
 * (strings) to file content (content_type).
 */
class ram_fs
{
    public:
        typedef std::vector<char> content_type;

    private:
        friend class ram_fs_initializer;
        typedef std::map<std::string, content_type> mss_type;
        static mss_type m_map;
        static std::recursive_mutex m_rlock;

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
};

//! Determines if the given file is a RAM-file.
bool is_ram_file(const std::string& file);

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
