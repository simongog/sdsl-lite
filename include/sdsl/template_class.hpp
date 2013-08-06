/*!\file template_class.hpp
   \brief template_class.hpp contains a template for our own sdsl class
   \author Your name
*/
#ifndef INCLUDED_SDSL_TEMPLATE_CLASS // replace TEMPLATE_CLASS by our class name
#define INCLUDED_SDSL_TEMPLATE_CLASS // in this include guard

#include "bit_vectors.hpp" // include your required sdsl header here
#include <iostream>        // then include other header, e.g. from the STL

//! Namespace for the succinct data structure library
namespace sdsl
{

//! Short description for what template_class is used here.
/*!
 *  Longer description of the class. If you implement a data structure from a
 *  paper, please mention this here in a separate references paragraph like this"
 *  \par References
 *       - Hugo Hacker: ``The hyper succinct data structure'', 2023.
 *
 *  Explain also the template parameters.
 *  \tparam bit_vector_type This template_class is dependent on one template parameter which
 *                         specifies the type of the bitvector which is used for a member variable.
 *                         This enables us to use the class with compressed and uncompressed
 *                         bitvector representation.
 */
template<class bit_vector_type=bit_vector>
class template_class
{
    public:
        typedef bit_vector::size_type size_type; // put public typedefs here
    private:
        size_type 		m_size; // member variables are prefixed with ``m_''
        bit_vector_type	m_data; // use descriptive names for the variables

    public:
        // Define a default constructor. Don't forget to initialize all members,
        // which are not initialized by default. E.g. m_size in this case.
        template_class():m_size(0) {}

        // User defined operator.
        explicit template_class(bit_vector b):m_size(b.size()) {
            bit_vector data(b.size()/2, 0);
            for (size_type i=0; i+1 < b.size(); i+=2) {
                data[i] = b[2*i] & b[2*i+1];
            }
            util::assign(m_data, data);
        }

        // copy constructor
        template_class(const template_class& t) {
            if (this != &t) { // if it is not the same object
                m_size = t.m_size;
                m_data = t.m_data;
            }
        }

        // You should implement the swap operator. It is used to swap two variables
        // of type template_class in constant time.
        void swap(const template_class& t) {
            std::swap(m_size, t.m_size); // use std::swap to swap variables of basic primitive types
            m_data.swap(t.m_data); // use the swap method for composite data types like sdsl or stl classes
        }

        // The serialize method writes the data structure into a stream ``out'' and returns
        // the number of written bytes.
        size_type serialize(std::ostream& out)const {
            size_type written_bytes = 0;
            // variables of basic primitive types can be written using write_member
            written_bytes += write_member(m_size, out);
            // other sdsl classes can be written by calling their serialize method
            written_bytes += m_data.serialize(out);
            return written_bytes;
        }

        // The load method reads the data structure from a stream ``in''. Note that
        // the members have to be read in the same order in which they were written.
        void load(std::istream& in) {
            // variables of basic primitive types can be read using read_member
            read_member(m_size, in);
            // other sdsl classes can be read by calling their load method
            m_data.load(in);
        }

        // user defined methods or operators
        bool operator[](size_type i)const {
            return m_data[i/2];
        }
};

}
#endif
