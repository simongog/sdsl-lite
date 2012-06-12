/*!\file structure_tree.hpp
   \brief structure_tree.hpp contains a helper class which can represent the memory structure of a class.
   \author Simon Gog
*/
#ifndef INCLUDED_SDSL_STRUCTURE_TREE
#define INCLUDED_SDSL_STRUCTURE_TREE

#include <vector>
#include <map>
#include <string>
#include <iostream>

using std::vector;
using std::map;
using std::string;

//! Namespace for the succinct data structure library
namespace sdsl
{

class structure_tree; // forward declaration

namespace util
{

template<typename T>
std::string to_string(const T&); // forward declaration

}

//! Class for a node of the structure tree
class structure_tree_node
{
    public:
        typedef map<string, string> tKeyValue;
        friend class structure_tree;

        structure_tree_node*&			 parent;
        vector<structure_tree_node*>&	 children;
        map<string, string>&			 key_values;
        //! Standard constructor
        structure_tree_node();
        //! User defined operator.
        /*! Adds the node as a child of node v.
         *  \param v The parent node of the created node.
         */
        explicit structure_tree_node(structure_tree_node* v);
        //! User defined operator.
        /*! Adds the node as a child of node v and sets two key values pairs
         *	\param v 			The parent node of the created node.
         *  \param name 		The name of the object.
         *  \param class_name   The class type of the object.
         */
        structure_tree_node(structure_tree_node* v, const string& name, const string& class_name);
        //! Copy constructor.
        structure_tree_node(const structure_tree_node& v);
        //! Destructor.
        ~structure_tree_node();
        structure_tree_node& operator=(const structure_tree_node& v);
        //! Swap operator.
        void swap(structure_tree_node& v);
        //! Add a key value pair.
        void add_key_value(const string& key, const string& value);
        template<class IntType>
        void add_size(IntType value);

    private:
        structure_tree_node*	 	 m_parent;
        vector<structure_tree_node*> m_children;
        map<string, string>			 m_key_values;
        void copy(const structure_tree_node& v);
};

class structure_tree
{
    public:
        static structure_tree_node* add_child(structure_tree_node* v, const string& name, const string& class_name);
        template<class IntType>
        static void add_size(structure_tree_node* v, IntType value);
        static structure_tree_node* parent(const structure_tree_node* v);
};


enum format_type {JSON_FORMAT, R_FORMAT};

template<format_type F>
void write_structure_tree(const structure_tree_node* v, std::ostream& out);


template<class IntType>
void structure_tree_node::add_size(IntType value)
{
    m_key_values["size"] = util::to_string(value);
}

template<class IntType>
void structure_tree::add_size(structure_tree_node* v, IntType value)
{
    if (NULL != v) {
        v->add_size(value);
    }
}

}
#endif
