#include <vector>
#include <stxxl/vector>

using namespace std;

int main()
{

    std::string asd = "21333213\t213213";

    istringstream ss(asd);

    uint var1, var2;

    ss >> var1 >> var2;

    std::cout << "var1:" << var1 << " var2: " << var2;

    return 0;
}
