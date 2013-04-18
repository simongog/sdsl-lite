#!/bin/bash

if [ $# -lt 1 ]; then
	echo "Usage: ${0} input tikz=0 divSpace=1 Time=seconds ylabelAdd=\" (in bytes))\""
	echo " input   : file where the output of mm-log.x is stored"
	echo " tikz = 0: output in mm-log.pdf"
	echo " tikz = 1: output in mm-log.tex"
	echo " divSpace: divide space by this number"
	echo " divTime : seconds | milliseconds"
	exit
fi

CUR_DIR=`pwd`

cd `dirname $1`
INPUT_DIR=`pwd`
INPUT=${INPUT_DIR}/`basename $1`

echo "INPUT=${INPUT}"
cd ${CUR_DIR}

OLD_DIR="$( cd "$( dirname "$0" )" && pwd )" # gets the directory where the script is located in
cd "${OLD_DIR}"

TIKZ="0"
if [ $# -gt 1 ]; then
	TIKZ=$2
fi
DIVSPACE=1
if [ $# -gt 2 ]; then
	DIVSPACE=$3
fi
TIME=seconds
if [ $# -gt 3 ]; then
	TIME=$4
fi
YLABELADD=" (in bytes)"
if [ $# -gt 4 ]; then
	YLABELADD=$5
fi

MM_DIR=`pwd`
echo "MM_DIR=${MM_DIR}"
R --vanilla --args "${MM_DIR}" "${INPUT}" "${TIKZ}" "${DIVSPACE}" "${TIME}" "${YLABELADD}" < mm-log.R
