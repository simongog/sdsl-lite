# find divsufsort library of Yuta Mori

# search in /usr/local/include and ~/include for divsufsort includes
find_path(divsufsort_INCLUDE_DIRS divsufsort.h /usr/local/include ~/include "$ENV{divsufsort_ROOT}")
# search in /usr/local/lib and ~/lib for divsufsort library 
find_library(divsufsort_LIBRARIES divsufsort /usr/local/lib ~/lib "$ENV{divsufsort_ROOT}")

set(divsufsort_FOUND TRUE)

if(NOT divsufsort_INCLUDE_DIRS)	
	set(divsufsort_FOUND FALSE)
endif(NOT divsufsort_INCLUDE_DIRS)	

if(NOT divsufsort_LIBRARIES)	
	set(divsufsort_FOUND FALSE)
else()	
	get_filename_component(divsufsort_LIBRARY_DIRS ${divsufsort_LIBRARIES} PATH)
	message("FOUND divsufsort library in ${divsufsort_LIBRARY_DIRS}")
endif(NOT divsufsort_LIBRARIES)	
