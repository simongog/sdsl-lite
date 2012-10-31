/* sdsl - succinct data structures library
    Copyright (C) 2008 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file testutils.hpp
 *  \brief testutils.hpp contains a "stopwatch" class for performance meassurement of program pieces.
 *  \author Simon Gog
 */
#ifndef INCLUDE_SDSL_TESTUTILS
#define INCLUDE_SDSL_TESTUTILS

#include "util.hpp"
#include "uintx_t.hpp"
#include <sys/time.h> // for struct timeval
#include <sys/resource.h> // for struct rusage
#include <iomanip>
#include <iostream>
#include <string>

namespace sdsl
{

//! A helper class to meassure the time consumption of program pieces.
/*! stop_watch is a stopwatch based on the commands getrusage and
 *  gettimeofday. Where getrusage is used to determine the user and system time
 *  and gettimeofday to determine the elapsed real time.
 */
class stop_watch
{
    private:
        rusage m_ruse1, m_ruse2;
        timeval m_timeOfDay1, m_timeOfDay2;
        static timeval m_first_t;
        static rusage m_first_r;
    public:

        stop_watch() : m_ruse1(), m_ruse2(), m_timeOfDay1(), m_timeOfDay2() {
            timeval t;
            t.tv_sec = 0; t.tv_usec = 0;
            m_ruse1.ru_utime = t; m_ruse1.ru_stime = t; // init m_ruse1
            m_ruse2.ru_utime = t; m_ruse2.ru_stime = t; // init m_ruse2
            m_timeOfDay1 = t; m_timeOfDay2 = t;
            if (m_first_t.tv_sec == 0) {
                gettimeofday(&m_first_t, 0);
            }
            if (m_first_r.ru_utime.tv_sec == 0 and m_first_r.ru_utime.tv_usec ==0) {
                getrusage(RUSAGE_SELF, &m_first_r);
            }
        }
        //! Start the stopwatch.
        /*! \sa stop
         */
        void start();

        //! Stop the stopwatch.
        /*! \sa start
         */
        void stop();

        //! Get the elapsed user time in milliseconds between start and stop.
        /*! \sa start, stop, get_real_time, get_sys_time
         */
        double get_user_time();

        //! Get the elapsed system time in milliseconds between start and stop.
        /*! \sa start, stop, get_real_time, get_user_time
         */
        double get_sys_time();

        //! Get the elapsed real time in milliseconds between start and stop.
        /*! \sa start, stop, get_sys_time, get_user_time
         */
        double get_real_time();

        //! Get the elapsed user time in milliseconds since the first construction of a stop_watch in the current process.
        /*! \sa get_user_time
         */
        uint64_t get_abs_user_time();

        //! Get the elapsed system time in milliseconds since the first construction of a stop_watch in the current process.
        /*! \sa get_sys_time
         */
        uint64_t get_abs_sys_time();

        //! Get the elapsed real time in milliseconds since the first construction of a stop_watch in the current process.
        /*! \sa get_real_time
         */
        uint64_t get_abs_real_time();

        uint64_t get_abs_page_faults();
};

//! Write stopwatch output in readable format
inline void write_R_output(std::string data_structure, std::string action,
                           std::string state="begin", uint64_t times=1, uint64_t check=0)
{
    if (util::verbose) {
        stop_watch _sw;
        _sw.stop();
        std::cout << data_structure << "\t" << action << "\t" << state << "\t"
                  << std::setw(9)<< times << "\t" << std::setw(9) << check << "\t"
                  << std::setw(9) << _sw.get_abs_real_time() << "\t "
                  << std::setw(9) << _sw.get_abs_user_time() << "\t"
                  << std::setw(9) << _sw.get_abs_sys_time() << std::endl;
    }
}

//! A helper class to get time information.
class clock
{
    public:
        static std::string get_time_string();
};



//! A helper class to handle files.
class file
{
    public:
        //! Read the file with the given file_name
        /*! \param file_name The file name of the text to read.
         *  \param c A char pointer which will point to the text that was read.
         *           New memory is allocated for the text. So free c if `read_text`
         *           was successful and c is not needed anymore.
         *  \param trunc Indicated if the file should be truncated.
         *  \param lim Maximal number of bytes which are read when trunc is true.
         *	\return len The number of read bits. If this is zero, now memory is
         *          allocated for c. And c equals NULL.
         *  \pre c has to be initialized to NULL.
         *  \post If len > 0  c[len]=0 and the memory for c was allocated with "new" else c=NULL.
         */
        static uint64_t read_text(const char* file_name, char*& c, bool trunc=0, uint64_t lim=0);

        static void write_text(const char* file_name, const char* c, uint64_t len);
};

} // end of namespace sdsl

#endif
