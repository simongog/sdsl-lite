/*! \file ram_fs.hpp
 * \brief ram_fs.hpp
 * \author Simon Gog
 */
#ifndef INCLUDED_SDSL_RAM_FS
#define INCLUDED_SDSL_RAM_FS

#include <string>
#include <map>

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
 * A RAM-file is represented as a (name, content)-pair.
 */
class ram_fs
{
        friend class ram_fs_initializer;
        typedef std::map<std::string, std::string> mss_type;
        static mss_type m_map;

    public:
        //! Default construct
        ram_fs();
        //! Store data under key `name`
        static void store(const std::string& name, const std::string& data);
        //! Get the content
        static const std::string& content(const std::string& name);
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
