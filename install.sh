#!/bin/bash
# This script builds all dependencies of sdsl
# and installs the library on a LINUX or Mac OS X system

CUR_DIR=`pwd`
SDSL_INSTALL_PREFIX=${HOME}
if [ $# -ge 1 ]; then
	SDSL_INSTALL_PREFIX=${1}
fi

# Get absolute path name of install directory
mkdir -p "${SDSL_INSTALL_PREFIX}" 2> /dev/null
cd "${SDSL_INSTALL_PREFIX}" > /dev/null 2>&1
if [ $? != 0 ] ; then
	echo "ERROR: directory '${SDSL_INSTALL_PREFIX}' does not exist nor could be created."
	echo "Please choose another directory."
	exit 1
else
	SDSL_INSTALL_PREFIX=`pwd -P`
fi

echo "Library will be installed in '${SDSL_INSTALL_PREFIX}'"

cd "${CUR_DIR}"
OLD_DIR="$( cd "$( dirname "$0" )" && pwd )" # gets the directory where the script is located in
cd "${OLD_DIR}"
OLD_DIR=`pwd`

# (1) Copy pre-commit hook


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

cd build # change into the build directory
if [ $? != 0 ]; then
	exit 1
fi
./clean.sh # clean-up build directory
if [ $? != 0 ]; then
	exit 1
fi

cmake -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_INSTALL_PREFIX="${SDSL_INSTALL_PREFIX}" .. # run cmake 
if [ $? != 0 ]; then
	echo "ERROR: CMake build failed."
	exit 1
fi
make sdsl # run make
if [ $? != 0 ]; then
	echo "ERROR: Build failed."
	exit 1
fi
echo "Removing old files"
echo "rm -rf '${SDSL_INSTALL_PREFIX}/include/sdsl/*'"
rm -rf "${SDSL_INSTALL_PREFIX}/include/sdsl/*"
if [ $? != 0 ]; then
	echo "WARNING: Could not remove old header files."
fi
echo "rm -f '${SDSL_INSTALL_PREFIX}/lib/libsdsl*'"
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
echo "The sdsl include files are located in '${SDSL_INSTALL_PREFIX}/include'."
echo "The library files are located in '${SDSL_INSTALL_PREFIX}/lib'."
echo " "
echo "Sample programs can be found in the examples-directory."
echo "A program 'example.cpp' can be compiled with the command: "
echo "g++ -std=c++11 -DNDEBUG -O3 [-msse4.2] \\"
echo "   -I${SDSL_INSTALL_PREFIX}/include -L${SDSL_INSTALL_PREFIX}/lib \\"
echo "   example.cpp -lsdsl -ldivsufsort -ldivsufsort64"
echo " "
echo "Tests in the test-directory"
echo "A cheat sheet in the extras/cheatsheet-directory."
echo "Have fun!"
