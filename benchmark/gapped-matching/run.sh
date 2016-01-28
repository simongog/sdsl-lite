#!/bin/bash

echo -e "\e[7m\e[1mBUILDING INDICES\e[0m"
make indices -j9
echo -e "\e[7m\e[1mEXTRACTING PHRASES\e[0m"
make inputs -j9 &> /dev/null
echo -e "\e[7m\e[1mRUN EXPERIMENTS\e[0m"
make timing
