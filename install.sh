#!/bin/bash
# This script downloads all dependencies of sdsl
# and installs the library on a LINUX or Mac OS system
#

# (1) Download libdivsufsort (excellent suffix array construction implementation by Yuta Mori) 
# version 2.0.1 from google code

OLD_DIR="$( cd "$( dirname "$0" )" && pwd )" # gets the directory where the script is located in
cd ${OLD_DIR}

source "scripts/download.sh"

LIBDIVSUFSORT_URL="http://libdivsufsort.googlecode.com/files/libdivsufsort-2.0.1.tar.gz"
LIBDIVSUFSORT_NAME=`basename ${LIBDIVSUFSORT_URL}`
if [[ ! -e ${LIBDIVSUFSORT_NAME} ]]; then
	echo "Download ${LIBDIVSUFSORT_NAME} from ${LIBDIVSUFSORT_URL}"
	download_from_url "${LIBDIVSUFSORT_URL}"
	if [[ $? != 0 ]] ; then  # if download_from_url returns not 0 
		echo "ERROR: Could not download the libdivsufsort library. Has the library moved to another place?"
		exit 1
	fi
fi

LIBDIVSUFSORT_DIR=${LIBDIVSUFSORT_NAME%.tar.gz}

if [[ ! -d "${LIBDIVSUFSORT_DIR}" ]]; then # if libdivsufsort is not unpacked
	echo "unpacking ${LIBDIVSUFSORT_NAME}"
	tar -xzvf ${LIBDIVSUFSORT_NAME} # 
	if [[ $? != 0 ]] ; then
		echo "ERROR: Could not unpack the libdivsufsort library."
		exit 1
	fi
fi

echo "LIBDIVSUFSORT_DIR=${LIBDIVSUFSORT_DIR}"
cd ${LIBDIVSUFSORT_DIR} # change in the libdivsufsort directory
echo ""
if [[ ! -d build ]]; then
	mkdir build # create a build directory
fi

cd build  # step into the build dir of libdivsufsort
cmake -DCMAKE_INSTALL_PREFIX=${HOME} -DBUILD_DIVSUFSORT64=ON \
	  -DBUILD_SHARED_LIBS=OFF -DBUILD_EXAMPLES=OFF .. # run cmake 
make # run make

if [[ ! -d "${HOME}/include"  ]]; then # if there exists not include dir 
	mkdir ${HOME}/include              # try to create one
fi
if [[ ! -d "${HOME}/lib" ]]; then     #if there exits no lib dir
	mkdir ${HOME}/lib
fi
make install  # install libdivsufsort

cd ../.. # go to the original directory
if [[ `pwd` != ${OLD_DIR} ]]; then
	echo "ERROR: we are not in the original dir ${OLD_DIR} now."
	exit 1
fi

echo "SUCCESS: Yuta Mori's libdivsufsort was installed successfully!"
rm -rf ${LIBDIVSUFSORT_DIR}
# (2) Download the google's test framework gtest 
#

GTEST_URL="http://googletest.googlecode.com/files/gtest-1.6.0.zip"
GTEST_NAME=`basename ${GTEST_URL}`

if [[ ! -e ${GTEST_NAME} ]]; then
	echo "Download ${GTEST_NAME} from ${GTEST_URL}"
	download_from_url "${GTEST_URL}"
	if [[ $? != 0 ]] ; then  # if curl returns not 0 
		echo "ERROR: Could not download the gtest library. Has the library moved to another place?"
		exit 1
	fi
fi
GTEST_DIR=${GTEST_NAME%.zip}
if [[ ! -d "${GTEST_DIR}" ]]; then # if gtest is not unpacked
	echo "unpacking ${GTEST_NAME}"
	unzip ${GTEST_NAME} # 
	if [[ $? != 0 ]] ; then
		echo "ERROR: Could not unpack the gtest library."
		exit 1
	fi
fi

cd ${GTEST_DIR} # change into the gtest directory
if [[ $? != 0 ]]; then
	exit 1
fi
if [[ ! -d build ]]; then
	mkdir build # create a build directory
fi
cd build  # step into the build dir of libdivsufsort
if [[ $? != 0 ]]; then
	exit 1
fi
cmake .. # run cmake
if [[ $? != 0 ]]; then
	exit 1
fi
make # run make
if [[ $? != 0 ]]; then
	exit 1
fi
cp libgtest.a libgtest_main.a ${HOME}/lib
cd ..
cp -r include/gtest ${HOME}/include
if [[ $? != 0 ]]; then
	exit 1
fi

cd ..

if [[ `pwd` != ${OLD_DIR} ]]; then
	echo "ERROR: we are not in the original dir ${OLD_DIR} now."
	exit 1
fi

echo "SUCCESS: gtest was installed successfully!"
rm -rf ${GTEST_DIR}

# (3) Finally, install sdsl
#

cd build # change into the build directory
if [[ $? != 0 ]]; then
	exit 1
fi

cmake -DCMAKE_INSTALL_PREFIX=${HOME} .. # run cmake 
make # run make
make install # install library

cd ..

if [[ `pwd` != ${OLD_DIR} ]]; then
	echo "ERROR: we are not in the original dir ${OLD_DIR} now."
	exit 1
fi

echo "SUCCESS: sdsl was installed successfully!"
echo "The sdsl include files are located in ${HOME}/include."
echo "The library files are located in ${HOME}/lib."
echo " "
echo "A program 'example.cpp' can be compiled with"
echo "g++ -DNDEBUG -O3 [-msse4.2] -I${HOME}/include -L${HOME}/lib example.cpp -lsdsl -ldivsufsort -ldivsufsort64"
echo " "
echo "Tests can be found in the test-directory."
echo "Have fun!"
