#include <sdsl/suffix_trees.hpp>
#include <iostream>

using namespace sdsl;
using namespace std;

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

int main(int argc, char** argv)
{
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " file" << endl;
        cout << " Creates a CST and CSA for a byte file and visualizes the memory utilization during construction." << endl;
        return 1;
    }

    if(0){
        memory_monitor::start();

        csa_sada<> csa;
        auto start = timer::now();
        construct(csa, argv[1], 1);
        auto stop = timer::now();
        cout << "construction csa time in seconds: " << duration_cast<seconds>(stop-start).count() << endl;

        memory_monitor::stop();
        std::ofstream csaofs("csa-construction_file.html");
        cout << "writing memory usage visualization to csa-construction.html\n";
        memory_monitor::write_memory_log<HTML_FORMAT>(csaofs);
        csaofs.close();
    }

    {
        //read file
        std::string text;
        {
            std::ostringstream ss;
            std::ifstream ifs(argv[1]);
            ss << ifs.rdbuf();
            text = ss.str();
        }

        memory_monitor::start();

        csa_sada<> csa;
        auto start = timer::now();
        /*
        {
            construct_im(csa, std::move(text), 1);
        }
        */
        construct(csa, argv[1], 1);
        auto stop = timer::now();
        cout << "construction csa time in seconds: " << duration_cast<seconds>(stop-start).count() << endl;

        memory_monitor::stop();
        std::ofstream csaofs("csa-construction_im.html");
        cout << "writing memory usage visualization to csa-construction.html\n";
        memory_monitor::write_memory_log<HTML_FORMAT>(csaofs);
        csaofs.close();
    }

}
