#!/bin/bash

mkdir collections || echo "found collections directory"
echo -e "\e[7m\e[1mCOMPILING\e[0m"
make execs
echo -e "\e[7m\e[1mBUILDING INDICES\e[0m"
make indices
echo -e "\e[7m\e[1mEXTRACTING PHRASES\e[0m"
make inputs -j9 &> /dev/null
echo -e "\e[7m\e[1mRUN EXPERIMENTS\e[0m"
make timing
