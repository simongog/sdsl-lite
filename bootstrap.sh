#! /bin/bash

username=d056848
################################get git
cd /local/${username}/
wget --no-check-certificate https://github.com/git/git/archive/master.zip
unzip master.zip
cd git-master
make configure ;# as yourself
mkdir -p /local/${username}/git-master/build
./configure --prefix=/local/d056848/git-master/build ;# as yourself
make all ;# as yourself
make install install-doc install-html;# as root
export PATH=/local/${username}/git-master/build/bin:$PATH

################################ get cmake
cd /local/${username}/
wget https://cmake.org/files/v3.5/cmake-3.5.2.tar.gz --no-check-certificate
tar -xzvf cmake-3.5.2.tar.gz 
cd cmake-3.5.2/
mkdir build
./bootstrap --prefix=/local/${username}/cmake-3.5.2/build && make && make install
export PATH=/local/${username}/cmake-3.5.2/build/bin:$PATH

################################ get gcc
cd /local/${username}/
wget ftp://ftp.mpi-sb.mpg.de/pub/gnu/mirror/gcc.gnu.org/pub/gcc/releases/gcc-4.9.3/gcc-4.9.3.tar.gz
cd gcc-4.9.3/
./contrib/download_prerequisites 
mkdir build
cd build
$PWD/../configure --prefix=$PWD --disable-multilib --enable-languages=c,c++,fortran,go
make -j 48
make install
export PATH=${PWD}/bin:$PATH


############################### install htop
cd /local/${username}/
wget http://hisham.hm/htop/releases/2.0.1/htop-2.0.1.tar.gz
tar -xzvf htop-2.0.1.tar.gz 
cd htop-2.0.1/
mkdir build
cd build
$PWD/../configure --prefix=$PWD
make
make install 
export PATH=${PWD}/bin:$PATH

################################ install bzip2 headers (for boost)
cd /local/${username}/
wget http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz
tar -xzvf bzip2-1.0.6.tar.gz 
cd bzip2-1.0.6/
mkdir build
make -f Makefile-libbz2_so
make install PREFIX=$PWD/build/
export PATH=${PWD}/build/bin:$PATH
export LD_LIBRARY_PATH=${PWD}/build/include:$PATH

################################ clone my sdsl-lite
cd /local/${username}/
GIT_SSL_NO_VERIFY=true git clone https://github.wdf.sap.corp/D056848/sdsl-lite.git
git config http.sslVerify "false"
git submodule update --init
cd build
rm -rf CMake* && cmake -D CMAKE_BUILD_TYPE=Release ..
make boost_external

echo -e "username=d056848
export PATH=/local/\${username}/git-master/build/bin:\$PATH
export PATH=/local/\${username}/cmake-3.5.2/build/bin:\$PATH
export PATH=/local/\${username}/gcc-4.9.3/build/bin:\$PATH
export PATH=/local/\${username}/htop-2.0.1/build/bin:\$PATH
export LIBRARY_PATH=/local/\${username}/bzip2-1.0.6/build/lib:\$LIBRARY_PATH
export PATH=/local/\${username}/bzip2-1.0.6/build/bin/:\$PATH
export CPLUS_INCLUDE_PATH=/local/\${username}/bzip2-1.0.6/build/include/:\$CPLUS_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=/local/\${username}/gcc-4.9.3/build/include/:\$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=/local/\${username}/gcc-4.9.3/build/lib:\$LIBRARY_PATH
export LIBRARY_PATH=/local/\${username}/gcc-4.9.3/build/lib64:\$LIBRARY_PATH
export LD_LIBRARY_PATH=/local/\${username}/gcc-4.9.3/build/lib:\$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/local/\${username}/gcc-4.9.3/build/lib64:\$LD_LIBRARY_PATH
"
>> ~/.bashrc
