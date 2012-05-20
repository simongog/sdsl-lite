#!/bin/bash

# Helper function to download files. It tries
# to use either wget or curl.
# Argument 1: URL of the download file
function download_from_url {
	wget_path=`which wget` # returns an empty path if wget is not found
	if [[ -z ${wget_path} ]] ; then # if empty path => wget not found
		curl_path=`which curl`
		if [[ -z ${curl_path} ]] ; then
			echo "ERROR: Neither wget nor curl is installed. Can not download required files."
		else # curl found => download with curl
			curl -O $1
		fi
	else # wget found => download with wget
		wget $1 
	fi	
}
