#!/bin/bash

CUR_DIR=`pwd`
OLD_DIR="$( cd "$( dirname "$0" )" && pwd )" # gets the directory where the script is located in
cd "${OLD_DIR}"

if [ $# -lt 1 ]; then
	echo "Usage: ${0} input tikz=0 divSpace=1 Time=seconds ylabelAdd=\" (in bytes))\""
	echo " input   : file where the output of mm-log.x is stored"
	echo " tikz = 0: output in mm-log.pdf"
	echo " tikz = 1: output in mm-log.tex"
	echo " divSpace: divide space by this number"
	echo " divTime : seconds|milliseconds"
	exit
fi

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

R --vanilla --args "../../benchmark/basic_functions.R" "$1" "${TIKZ}" "${DIVSPACE}" "${TIME}" "${YLABELADD}" < mm-log.R


