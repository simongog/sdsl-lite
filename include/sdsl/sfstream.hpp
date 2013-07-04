/*!\file sfstream.hpp
   \brief sfstream.hpp contains a two stream class which can be used to read/write from/to files or strings.
   \author Simon Gog
*/
#ifndef INCLUDED_SDSL_SFSTREAM
#define INCLUDED_SDSL_SFSTREAM

#include <fstream>
#include <sstream>
#include <string>

namespace sdsl
{

class osfstream : public std::ostream
{
    private:
        std::streambuf* m_streambuf;
        std::string      m_file;
        bool             m_use_ram;
        bool             m_closed; // indicates if the buffer is closed
    public:
        typedef void* voidptr;
        //! Standard constructor.
        osfstream();
        //! Constructor taking a file name and open mode.
        osfstream(const std::string& file, std::ios_base::openmode mode = std::ios_base::out);
        //! Open the stream.
        std::streambuf*
        open(const std::string& file, std::ios_base::openmode mode = std::ios_base::out);
        //! Is the stream close?
        bool is_open();
        //! Close the stream.
        void close();
        //! Standard destructor
        ~osfstream();
        //! Cast to void*
        operator  voidptr() const;
};


class isfstream : public std::istream
{
    private:
        std::streambuf* m_streambuf;
        std::string      m_file;
        bool             m_use_ram;
        bool             m_closed; // indicates if the buffer is closed
    public:
        typedef void* voidptr;
        //! Standard constructor.
        isfstream();
        //! Constructor taking a file name and open mode.
        isfstream(const std::string& file, std::ios_base::openmode mode = std::ios_base::in);
        //! Open the stream.
        std::streambuf*
        open(const std::string& file, std::ios_base::openmode mode = std::ios_base::in);
        //! Is the stream close?
        bool is_open();
        //! Close the stream.
        void close();
        //! Standard destructor
        ~isfstream();
        //! Cast to void*
        operator  voidptr() const;
};

} // end namespace

#endif
