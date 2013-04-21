#include <iostream>
#include <string>
#include <sstream>

using namespace std;

const int BUF_SIZE=4096;
char line[BUF_SIZE];

int main()
{
    string token;
    string pat="candidates=";
    cin.getline(line, BUF_SIZE);
    stringstream ss(line);
    while (ss >> token) {
        if (token.substr(0, pat.size()) == pat) {
            cout << token.substr(pat.size());
            return 0;
        }
    }
    cout << "-1" << endl;
}
