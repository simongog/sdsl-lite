echo
echo -e '\033[1m REBUILD \033[0m'
echo
g++ -std=c++11 -O3 -ffast-math -funroll-loops -DNDEBUG -I ../include -I ../build/external/libdivsufsort-2.0.1/include -L ~/lib scratch.cpp -o scratch -lsdsl -ldivsufsort -ldivsufsort64
./scratch input.txt