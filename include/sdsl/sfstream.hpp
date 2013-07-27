/*!\file sfstream.hpp
   \brief sfstream.hpp contains a two stream class which can be used to read/write from/to files or strings.
   \author Simon Gog
*/
#ifndef INCLUDED_SDSL_SFSTREAM
#define INCLUDED_SDSL_SFSTREAM

#include <fstream>
#include <sstream>
#include <string>
#include "sdsl/ram_fs.hpp"
#include "sdsl/ram_filebuf.hpp"

namespace sdsl
{

class osfstream : public std::ostream
{
    public:
        typedef std::streambuf* buf_ptr_type;
    private:
        buf_ptr_type m_streambuf = nullptr;
        std::string  m_file      = "";
    public:
        typedef void* voidptr;
        //! Standard constructor.
        osfstream();
        //! Constructor taking a file name and open mode.
        osfstream(const std::string& file, std::ios_base::openmode mode = std::ios_base::out);
        //! Open the stream.
        buf_ptr_type
        open(const std::string& file, std::ios_base::openmode mode = std::ios_base::out);
        //! Is the stream close?
        bool is_open();
        //! Close the stream.
        void close();
        //! Standard destructor
        ~osfstream();
        //! Cast to void*
        operator  voidptr() const;

        osfstream& seekp(pos_type pos);
        osfstream& seekp(off_type off, ios_base::seekdir way);
        std::streampos tellp();
};


class isfstream : public std::istream
{
        typedef std::streambuf* buf_ptr_type;
    private:
        buf_ptr_type m_streambuf = nullptr;
        std::string  m_file      = "";
    public:
        typedef void* voidptr;
        //! Standard constructor.
        isfstream();
        //! Constructor taking a file name and open mode.
        isfstream(const std::string& file, std::ios_base::openmode mode = std::ios_base::in);
        //! Open the stream.
        buf_ptr_type
        open(const std::string& file, std::ios_base::openmode mode = std::ios_base::in);
        //! Is the stream close?
        bool is_open();
        //! Close the stream.
        void close();
        //! Standard destructor
        ~isfstream();
        //! Cast to void*
        operator  voidptr() const;

        isfstream& seekg(pos_type pos);
        isfstream& seekg(off_type off, ios_base::seekdir way);
        std::streampos tellg();
};

} // end namespace

#endif
