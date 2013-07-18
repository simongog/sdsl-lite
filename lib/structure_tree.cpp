#include "sdsl/structure_tree.hpp"
#include "sdsl/util.hpp"

namespace sdsl
{

void output_tab(std::ostream& out,size_t level)
{
    for (size_t i=0; i<level; i++) out << "\t";
}

template<>
void write_structure_tree<JSON_FORMAT>(const structure_tree_node* v, std::ostream& out, size_t level)
{
    if (v) {
        output_tab(out,level); out << "{" << std::endl;
        output_tab(out,level+1); out << "\"class_name\":" << "\"" << v->type << "\"," << std::endl;
        output_tab(out,level+1); out << "\"name\":" << "\"" << v->name << "\"," << std::endl;
        output_tab(out,level+1); out << "\"size\":" << "\"" << v->size << "\"";

        if (v->children.size()) {
            out << "," << std::endl; // terminate the size tag from before if there are children
            output_tab(out,level+1); out << "\"children\":[" << std::endl;
            size_t written_child_elements = 0;
            for (const auto& child : v->children) {
                if (written_child_elements++ > 0) {
                    out << "," << std::endl;
                }
                write_structure_tree<JSON_FORMAT>(child.second.get(), out,level+2);
            }
            out << std::endl;
            output_tab(out,level+1); out << "]" << std::endl;
        } else {
            out << std::endl;
        }
        output_tab(out,level); out << "}";
    }
}

template<>
void write_structure_tree<R_FORMAT>(const structure_tree_node* v, std::ostream& out, size_t level)
{
    /*    if (NULL == v or (v->children.size()==0 and v->key_values.size()==0)) {
            return;
        }
        typedef structure_tree_node::tKeyValue::const_iterator const_iterator;
        size_t written_elements = 0;
        out << "list("; // begin R list
        for (const_iterator it = v->key_values.begin(); it != v->key_values.end(); ++it) {
            if (written_elements++ > 0) {
                out << ",";
            }
            out << it->first << " = " << it->second;
        }
        if (v->children.size() > 0) {
            if (written_elements++ > 0) {
                out << ",";
            }
            out << "list("; // open children
            size_t written_child_elements = 0;
            for (size_t i = 0; i < v->children.size(); ++i) {
                if (written_child_elements++ > 0) {
                    out << ",";
                }
                write_structure_tree<JSON_FORMAT>(v->children[i], out);
            }
            out << ")";// close children
        }
        out << ")"; // end R list*/
}

} // namespace end
