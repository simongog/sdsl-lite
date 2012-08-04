#!/bin/bash
# This script builds all dependencies of sdsl
# and installs the library on a LINUX or Mac OS X system

SDSL_INSTALL_PREFIX=${HOME}
if [[ $# -ge 1 ]]; then
	SDSL_INSTALL_PREFIX=${1}
	echo "Library will be installed in ${1}"
fi


sufsort_dir="libdivsufsort-2.0.1"
gtest_dir="gtest-1.6.0"

# (1) Install libdivsufsort-2.0.1
OLD_DIR="$( cd "$( dirname "$0" )" && pwd )" # gets the directory where the script is located in
cd "${OLD_DIR}"

cd "${sufsort_dir}/build"  # step into the build dir of libdivsufsort

cmake -DCMAKE_INSTALL_PREFIX=${SDSL_INSTALL_PREFIX} -DBUILD_DIVSUFSORT64=ON \
	  -DBUILD_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE="Release" -DBUILD_EXAMPLES=OFF .. # run cmake 
make # run make

if [[ ! -d "${SDSL_INSTALL_PREFIX}/include"  ]]; then # if there exists not include dir 
	mkdir "${SDSL_INSTALL_PREFIX}/include"              # try to create one
fi
if [[ ! -d "${SDSL_INSTALL_PREFIX}/lib" ]]; then     #if there exits no lib dir
	mkdir "${SDSL_INSTALL_PREFIX}/lib"
fi
make install  # install libdivsufsort

cd ../.. # go to the original directory
if [[ `pwd` != ${OLD_DIR} ]]; then
	echo "ERROR: we are not in the original dir ${OLD_DIR} now."
	exit 1
fi
echo "SUCCESS: Yuta Mori's libdivsufsort was installed successfully!"
rm -rf "${sufsort_dir}/build/*"
# (2) Install gtest

cd "${gtest_dir}/build"  # step into the build dir of gtest

# (3) Finally, install sdsl
cmake -DCMAKE_INSTALL_PREFIX=${SDSL_INSTALL_PREFIX} \
	  -DBUILD_SHARED_LIBS=OFF -Dgtest_disable_pthreads=ON .. 
make # run make
cp ./libgtest* "${SDSL_INSTALL_PREFIX}/lib" # copy generated libraries to the destination
cp -r ../include/gtest "${SDSL_INSTALL_PREFIX}/include"

cd ../.. # go to the original directory
if [[ `pwd` != ${OLD_DIR} ]]; then
	echo "ERROR: we are not in the original dir ${OLD_DIR} now."
	exit 1
fi
echo "SUCCESS: gtest was installed successfully!"
#rm -rf "${sufsort_dir}/build/*"


cd build # change into the build directory
if [[ $? != 0 ]]; then
	exit 1
fi

cmake -DCMAKE_INSTALL_PREFIX=${SDSL_INSTALL_PREFIX} .. # run cmake 
make # run make
make install # install library

cd ..

if [[ `pwd` != ${OLD_DIR} ]]; then
	echo "ERROR: we are not in the original dir ${OLD_DIR} now."
	exit 1
fi

echo "SUCCESS: sdsl was installed successfully!"
echo "The sdsl include files are located in ${SDSL_INSTALL_PREFIX}/include."
echo "The library files are located in ${SDSL_INSTALL_PREFIX}/lib."
echo " "
echo "Sample programs can be found in the examples-directory."
echo "Test cases in the test-directory"
echo "A program 'example.cpp' can be compiled with the command: "
echo "g++ -DNDEBUG -O3 [-msse4.2] \\"
echo "   -I${SDSL_INSTALL_PREFIX}/include -L${SDSL_INSTALL_PREFIX}/lib \\"
echo "   example.cpp -lsdsl -ldivsufsort -ldivsufsort64"
echo " "
echo "Tests can be found in the test-directory."
echo "Have fun!"
