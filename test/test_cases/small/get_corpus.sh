#! /bin/bash

OLD_DIR=`pwd`
CUR_DIR="$( cd "$( dirname "$0" )" && pwd )" # gets the directory where the script is located in
cd ${CUR_DIR}

# Download test cases
GUTENBERG_MIRROR="ftp://www.artfiles.org/gutenberg.org"
FAUST_URL="${GUTENBERG_MIRROR}/2/2/2/2229/2229-8.txt"
FAUST_BASENAME=`basename ${FAUST_URL}`
FAUST_NAME='faust.txt'

ZARATHUSTRA_URL="${GUTENBERG_MIRROR}/7/2/0/7205/7205-8.txt"
ZARATHUSTRA_BASENAME=`basename ${ZARATHUSTRA_URL}`
ZARATHUSTRA_NAME='zarathustra.txt'

SIMONS_HOMEPAGE="http://people.eng.unimelb.edu.au/sgog/test_cases"
INT_TEST_CASES='keeper.int moby.int'

# Include download function
source "${CUR_DIR}/../../../scripts/download.sh"

if [[ ! -e ${FAUST_NAME} ]]; then
	download_from_url "${FAUST_URL}"
	mv ${FAUST_BASENAME} ${FAUST_NAME}
fi

if [[ ! -e ${ZARATHUSTRA_NAME} ]]; then
	download_from_url "${ZARATHUSTRA_URL}"
	mv ${ZARATHUSTRA_BASENAME} ${ZARATHUSTRA_NAME}
fi

for int_file in ${INT_TEST_CASES}; do
	if [[ ! -e int_file ]]; then
	download_from_url "${SIMONS_HOMEPAGE}/${int_file}"
	fi
done

cd ${OLD_DIR} 
