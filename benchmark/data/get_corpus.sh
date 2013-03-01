#! /bin/bash 

# Download and extract the files of the Pizza & Chili Corpus
# of P. Ferragina and G. Navarro
#
# Author: Simon Gog (simon.gog@unimelb.edu.au)

baseurl='http://pizzachili.dcc.uchile.cl/texts/'
files='nlang/english code/sources protein/proteins xml/dblp.xml dna/dna'
sizes='200'

my_dir="$( cd "$( dirname "$0" )" && pwd )" # gets the directory where the script is located in
cd "${my_dir}"

for file in ${files}; do
	for size in ${sizes}; do 
		file_name=${baseurl}${file}.${size}MB.gz
		file_name1=${baseurl}${file}.${size}MB
		basename=`basename ${file_name}`
		basename1=`basename ${file_name1}`
		echo ${file_name}
		echo ${basename}
# download the file if it does not exists on the disk
		if [[ (! -e ${basename1}) && (! -e ${basename}) ]]; then 
			curl -O ${file_name}
		fi	
# extract the file		
		if [ -e ${basename} ]; then
			gunzip ${basename}
		fi	
	done	
done
