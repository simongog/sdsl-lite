#include <sdsl/suffix_trees.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

std::string format("%3I%3S %3s %3P %3p %3L %3B  %T");
std::string header("  i SA ISA PSI  LF LCP BWT  TEXT");
static const std::streamsize BUF_SIZE=4096;
char text[BUF_SIZE];


typedef csa_wt<wt_int<>, 4, 4, sa_order_sa_sampling<>, int_vector<>, int_alphabet<> > csa_int_t;
typedef cst_sct3<>          cst_byte_t;
typedef cst_sct3<csa_int_t> cst_int_t;

void print_usage(const char* command)
{
    cout << "A pretty printer for suffix array/tree information." << endl;
    cout << "Usage: " << command << " X [format] [header]" << endl;
    cout << " X    : Input is interpreted dependend on X." << endl;
    cout << "        As byte sequence for X=1. " << endl;
    cout << "        As sequence of decimal numbers for X=d."<< endl;
    cout <<" format: Format string. Default=`" << format << "`." << endl;
    cout <<" header: Header string. Default=`" << header <<"`." << endl;
}

int main(int argc, char* argv[])
{
    if (argc < 2 or !('1' == argv[1][0] or 'd'==argv[1][0])) {
        print_usage(argv[0]);
        return 1;
    }
    if (argc > 2) {
        format = argv[2];
    }
    if (argc > 3) {
        header = argv[3];
    }
    while (cin.getline(text, BUF_SIZE)) {
        cout << header  << endl;
        if ('1' == argv[1][0]) {
            cst_byte_t cst;
            construct_im(cst, (const char*)text, 1);
            csXprintf(cout, format, cst, '$');
        } else if ('d' == argv[1][0]) {
            cst_int_t cst;
            construct_im(cst, (const char*)text, 'd');
            csXprintf(cout, format, cst, '0');
        }
        cout << endl;
    }
}
