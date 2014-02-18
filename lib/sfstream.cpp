#include "sdsl/sfstream.hpp"
#include "sdsl/util.hpp"
#include <iostream>

namespace sdsl
{

//  IMPLEMENTATION OF OSFSTREAM

osfstream::osfstream() : std::ostream(nullptr)
{
    this->init(m_streambuf);
}

osfstream::osfstream(const std::string& file, std::ios_base::openmode mode) : std::ostream(nullptr)
{
    this->init(m_streambuf);
    open(file, mode);
}

osfstream::buf_ptr_type
osfstream::open(const std::string& file, std::ios_base::openmode mode)
{
    delete m_streambuf;
    m_streambuf = nullptr;
    m_file = file;
    std::streambuf* success = nullptr;
    if (is_ram_file(file)) {
        m_streambuf = new ram_filebuf();
        success = ((ram_filebuf*)m_streambuf)->open(m_file, mode);
    } else {
        m_streambuf = new std::filebuf();
        success = ((std::filebuf*)m_streambuf)->open(m_file, mode);
    }
    if (success) {
        this->clear();
    } else {
        this->setstate(std::ios_base::failbit);
        delete m_streambuf;
        m_streambuf = nullptr;
    }
    this->rdbuf(m_streambuf);
    return m_streambuf;
}

bool
osfstream::is_open()
{
    if (nullptr == m_streambuf)
        return false;
    if (is_ram_file(m_file)) {
        return ((ram_filebuf*)m_streambuf)->is_open();
    } else {
        return ((std::filebuf*)m_streambuf)->is_open();
    }
}

void
osfstream::close()
{
    bool fail = false;
    if (nullptr == m_streambuf) {
        fail = true;
    } else {
        if (is_ram_file(m_file)) {
            fail = !((ram_filebuf*)m_streambuf)->close();
        } else {
            fail = !((std::filebuf*)m_streambuf)->close();
        }
    }
    if (fail) this->setstate(std::ios::failbit);
}

osfstream::~osfstream()
{
    delete m_streambuf; // streambuf closes the file on destruction
}

osfstream::operator voidptr()const
{
    return m_streambuf;
}

osfstream&
osfstream::seekp(pos_type pos)
{
    ios_base::iostate err = std::ios_base::iostate(std::ios_base::goodbit);
    try {
        if (!this->fail()) {
            pos_type p = 0;
            if (is_ram_file(m_file)) {
                p = ((ram_filebuf*)m_streambuf)->pubseekpos(pos, std::ios_base::out);
            } else {
                p = ((std::filebuf*)m_streambuf)->pubseekpos(pos, std::ios_base::out);
            }
            if (p == pos_type(off_type(-1))) {
                err |= ios_base::failbit;
                this->setstate(err);
            }
        }
    } catch (...) {
        if (err) {
            this->setstate(err);
        }
    }
    return *this;
}


osfstream&
osfstream::seekp(off_type off, std::ios_base::seekdir way)
{
    ios_base::iostate err = std::ios_base::iostate(ios_base::goodbit);
    try {
        if (!this->fail()) {
            pos_type p = 0;
            if (is_ram_file(m_file)) {
                p = ((ram_filebuf*)m_streambuf)->pubseekoff(off, way, std::ios_base::out);

            } else {
                p = ((std::filebuf*)m_streambuf)->pubseekoff(off, way, std::ios_base::out);
            }
            if (p == pos_type(off_type(-1))) {
                err |= ios_base::failbit;
                this->setstate(err);
            }
        }
    } catch (...) {
        if (err) {
            this->setstate(err);
        }
    }
    return *this;
}



//  IMPLEMENTATION OF ISFSTREAM

isfstream::isfstream() : std::istream(nullptr)
{
    this->init(m_streambuf);
}

isfstream::isfstream(const std::string& file, std::ios_base::openmode mode) : std::istream(nullptr)
{
    this->init(m_streambuf);
    open(file, mode);
}

isfstream::buf_ptr_type
isfstream::open(const std::string& file, std::ios_base::openmode mode)
{
    delete m_streambuf;
    m_streambuf = nullptr;
    m_file = file;
    std::streambuf* success = nullptr;
    if (is_ram_file(file)) {
        m_streambuf = new ram_filebuf();
        success = ((ram_filebuf*)m_streambuf)->open(m_file, mode);
    } else {
        m_streambuf = new std::filebuf();
        success = ((std::filebuf*)m_streambuf)->open(m_file, mode);
    }
    if (success) {
        this->clear();
    } else {
        this->setstate(std::ios_base::failbit);
        delete m_streambuf;
        m_streambuf = nullptr;
    }
    this->rdbuf(m_streambuf);
    return m_streambuf;
}

bool
isfstream::is_open()
{
    if (nullptr == m_streambuf)
        return false;
    if (is_ram_file(m_file)) {
        return ((ram_filebuf*)m_streambuf)->is_open();
    } else {
        return ((std::filebuf*)m_streambuf)->is_open();
    }
}

void
isfstream::close()
{
    bool fail = false;
    if (nullptr == m_streambuf) {
        fail = true;
    } else {
        if (is_ram_file(m_file)) {
            fail = !((ram_filebuf*)m_streambuf)->close();
        } else {
            fail = !((std::filebuf*)m_streambuf)->close();
        }
    }
    if (fail) this->setstate(std::ios::failbit);
}

isfstream&
isfstream::seekg(pos_type pos)
{
    ios_base::iostate err = std::ios_base::iostate(std::ios_base::goodbit);
    try {
        if (!this->fail()) {
            pos_type p = 0;
            if (is_ram_file(m_file)) {
                p = ((ram_filebuf*)m_streambuf)->pubseekpos(pos, std::ios_base::in);

            } else {
                p = ((std::filebuf*)m_streambuf)->pubseekpos(pos, std::ios_base::in);
            }
            if (p == pos_type(off_type(-1))) {
                err |= ios_base::failbit;
            }
        }
    } catch (...) {
        if (err) {
            this->setstate(err);
        }
    }
    return *this;
}


isfstream&
isfstream::seekg(off_type off, std::ios_base::seekdir way)
{
    ios_base::iostate err = std::ios_base::iostate(ios_base::goodbit);
    try {
        if (!this->fail()) {
            pos_type p = 0;
            if (is_ram_file(m_file)) {
                p = ((ram_filebuf*)m_streambuf)->pubseekoff(off, way, std::ios_base::in);

            } else {
                p = ((std::filebuf*)m_streambuf)->pubseekoff(off, way, std::ios_base::in);
            }
            if (p == pos_type(off_type(-1))) {
                err |= ios_base::failbit;
            }
        }
    } catch (...) {
        if (err) {
            this->setstate(err);
        }
    }
    return *this;
}

std::streampos
isfstream::tellg()
{
    ios_base::iostate err = std::ios_base::iostate(ios_base::goodbit);
    pos_type p = pos_type(off_type(-1));
    try {
        if (!this->fail()) {
            if (is_ram_file(m_file)) {
                p = ((ram_filebuf*)m_streambuf)->pubseekoff(0, std::ios_base::cur);

            } else {
                p = ((std::filebuf*)m_streambuf)->pubseekoff(0, std::ios_base::cur);
            }
            if (p == pos_type(off_type(-1))) {
                err |= ios_base::failbit;
            }
        }
    } catch (...) {
        if (err) {
            this->setstate(err);
        }
    }
    return p;
}

isfstream::~isfstream()
{
    delete m_streambuf;
}

isfstream::operator voidptr()const
{
    return m_streambuf; // streambuf closes the file on destruction
}

}// end namespace sdsl
