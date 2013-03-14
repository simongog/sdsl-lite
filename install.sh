#!/bin/bash
# This script builds all dependencies of sdsl
# and installs the library on a LINUX or Mac OS X system

SDSL_INSTALL_PREFIX=${HOME}
if [ $# -ge 1 ]; then
	SDSL_INSTALL_PREFIX=${1}
	echo "Library will be installed in ${1}"
fi

# (1) Copy pre-commit hook

OLD_DIR="$( cd "$( dirname "$0" )" && pwd )" # gets the directory where the script is located in
cd "${OLD_DIR}"

if [ -d ".git/hooks" ]; then
	echo "Copy pre-commit into .git/hooks"
	cp extras/pre-commit .git/hooks/
	if [ $? != 0 ]; then
		echo "WARNING: could not copy pre-commit script into .git/hooks"
	fi
	chmod u+x .git/hooks/pre-commit
	if [ $? != 0 ]; then
		echo "WARNING: could not make pre-commit script executable"
	fi
else
	echo "WARNING: .git/hooks directory does not exists."
	echo "         The pre-commit hook is not installed."
fi

# (2) Install divsufsort, gtest, and sdsl

cd build # change into the build directory
if [ $? != 0 ]; then
	exit 1
fi
rm -f CMakeCache.txt

cmake -DCMAKE_INSTALL_PREFIX=${SDSL_INSTALL_PREFIX} .. # run cmake 
if [ $? != 0 ]; then
	echo "ERROR: CMake build failed."
	exit 1
fi
make # run make
if [ $? != 0 ]; then
	echo "ERROR: Build failed."
	exit 1
fi
echo "Removing old files"
echo "rm -rf ${SDSL_INSTALL_PREFIX}/include/sdsl/*"
rm -rf "${SDSL_INSTALL_PREFIX}/include/sdsl/*"
if [ $? != 0 ]; then
	echo "WARNING: Could not remove old header files."
fi
echo "rm -f ${SDSL_INSTALL_PREFIX}/lib/libsdsl*"
rm -f "${SDSL_INSTALL_PREFIX}/lib/libsdsl*"
if [ $? != 0 ]; then
	echo "WARNING: Could not remove old library file."
fi
make install # install library
if [ $? != 0 ]; then
	echo "ERROR: Installation failed."
	exit 1
fi

cd ..

if [ "`pwd`" != "${OLD_DIR}" ]; then
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
