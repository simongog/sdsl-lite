#include <iostream>
#include <map>
#include <string>
#include <cstdint>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

const size_t buf_size=1024*128;
char   buffer[buf_size];

int main(int argc, char* argv[])
{
    if (argc < 3) {
        cout << "Usage: " << argv[0] << " pattern-file .dic-file" << endl;
        return 1;
    }
    ifstream dic_in(argv[2]);
    if (!dic_in) {
        cout << "Could not open dictionary file '" << argv[2] << "'" << endl;
        return 1;
    }
    ifstream pat_in(argv[1]);
    if (!pat_in) {
        cout << "Could not open pattern file '" << argv[1] << "'" << endl;
        return 1;
    }
    std::map<uint64_t, string> i2w;
    string word;
    while (dic_in >> word) {
        uint64_t nr, occ;
        dic_in >> nr >> occ;
        i2w[nr] = word;
    }
    while (pat_in.getline(buffer, buf_size)) {
        stringstream ss;
        ss << string(buffer);
        uint64_t w;
        while (ss >> w) {
            cout << i2w[w] << " ";
        }
        cout << "\n";
    }
}
