#!/bin/bash

CUR_DIR=`pwd`
cd "${CUR_DIR}"
OLD_DIR="$( cd "$( dirname "$0" )" && pwd )" # gets the directory where the script is located in
cd "${OLD_DIR}"
OLD_DIR=`pwd`

lcov --directory .. --zerocounters

cd ..
make clean
make CODE_COVER="-g -fprofile-arcs -ftest-coverage -lgcov" $1
MAKERET=$?
cd ${OLD_DIR}

if [ ${MAKERET} != 0 ]; then
	echo "ERROR: tests failed."
	exit 1
fi

INC_STR=`grep 'INC_DIR = ' ../../Make.helper`
INC_DIR=${INC_STR#INC_DIR = }

lcov -d .. -d ${INC_DIR}/sdsl --no-recursion --no-external --capture --output-file sdsl-lite.info
genhtml -o output --legend sdsl-lite.info
