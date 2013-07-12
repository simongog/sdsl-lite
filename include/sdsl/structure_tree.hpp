/*!\file structure_tree.hpp
   \brief structure_tree.hpp contains a helper class which can represent the memory structure of a class.
   \author Simon Gog
*/
#ifndef INCLUDED_SDSL_STRUCTURE_TREE
#define INCLUDED_SDSL_STRUCTURE_TREE

#include "uintx_t.hpp"
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <sstream>

using std::vector;
using std::map;

//! Namespace for the succinct data structure library
namespace sdsl
{

static const char STRUCTURE_TREE_SIZE_KEY[] = "size";
static const char STRUCTURE_TREE_NAME_KEY[] = "name";
static const char STRUCTURE_TREE_CLASS_NAME_KEY[] = "class_name";

class structure_tree; // forward declaration

//! Class for a node of the structure tree
class structure_tree_node
{
    public:
        typedef map<std::string, std::string> tKeyValue;
        friend class structure_tree;
    private:
        structure_tree_node*	 		m_parent;
        vector<structure_tree_node*> 	m_children;
        map<std::string, std::string> 	m_key_values;
        void copy(const structure_tree_node& v);
        void delete_children();
        bool equal_key_value_pairs(const tKeyValue& kv1, const tKeyValue& kv2)const;
    public:

        structure_tree_node*&			 parent;
        vector<structure_tree_node*>&	 children;
        map<std::string, std::string>&	 key_values;
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
        structure_tree_node(structure_tree_node* v, const std::string& name, const std::string& class_name);
        //! Copy constructor.
        structure_tree_node(const structure_tree_node& v);
        //! Destructor.
        ~structure_tree_node();
        structure_tree_node& operator=(const structure_tree_node& v);
        //! Swap operator.
        void swap(structure_tree_node& v);
        //! Add a key value pair.
        void add_key_value(const std::string& key, const std::string& value);
        void add_size(uint64_t value);

        //! Adds the size values of structure_tree_node v recursively to that of the current node.
        /*! \return True, if the recursive structure of both nodes is exactly the same and
         *          the size values were integers.
         */
        bool merge(const structure_tree_node& v);

    private:
        bool equal_structure(const structure_tree_node& v)const;

        void rec_merge(const structure_tree_node& v);
};

class structure_tree
{
    public:
        static structure_tree_node* add_child(structure_tree_node* v, const std::string& name, const std::string& class_name);
        static void add_size(structure_tree_node* v, uint64_t value);
        static structure_tree_node* parent(const structure_tree_node* v);
        static bool merge_children(structure_tree_node* v);
};


enum format_type {JSON_FORMAT, R_FORMAT};

template<format_type F>
void write_structure_tree(const structure_tree_node* v, std::ostream& out);


}
#endif
