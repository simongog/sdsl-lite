# find the 64 bit version divsufsort64	

find_path(divsufsort64_INCLUDE_DIRS divsufsort64.h /usr/local/include ~/include "$ENV{divsufsort64_ROOT}")
find_library(divsufsort64_LIBRARIES divsufsort64 /usr/local/lib ~/lib "$ENV{divsufsort64_ROOT}")

set(divsufsort64_FOUND TRUE)

if(NOT divsufsort64_INCLUDE_DIRS)	
	set(divsufsort64_FOUND FALSE)
endif(NOT divsufsort64_INCLUDE_DIRS)	

if(NOT divsufsort64_LIBRARIES)	
	set(divsufsort64_FOUND FALSE)
else()	
	get_filename_component(divsufsort64_LIBRARY_DIRS ${divsufsort_LIBRARIES} PATH)
	message("FOUND divsufsort library in ${divsufsort64_LIBRARY_DIRS}")

endif(NOT divsufsort64_LIBRARIES)	
