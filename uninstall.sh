#!/bin/bash
# This script removes all installed  sdsl files
# on a LINUX or Mac OS X system

CUR_DIR=`pwd`
SDSL_INSTALL_PREFIX=${HOME}
if [ $# -ge 1 ]; then
	SDSL_INSTALL_PREFIX=${1}
fi

echo "Library files will be removed from"
echo "'${SDSL_INSTALL_PREFIX}/lib' and"
echo "'${SDSL_INSTALL_PREFIX}/include'"


cd "${CUR_DIR}"
OLD_DIR="$( cd "$( dirname "$0" )" && pwd )" # gets the directory where the script is located in
cd "${OLD_DIR}"
OLD_DIR=`pwd`

cd build # change into the build directory
if [ $? != 0 ]; then
	exit 1
fi

make uninstall

if [ $? != 0 ]; then
    exit 1
fi

./clean.sh # clean-up build directory
if [ $? != 0 ]; then
	exit 1
fi

echo "Installed sdsl files were removed."
