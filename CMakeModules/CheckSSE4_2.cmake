# Check if the CPU provides fast operations
# for popcount, leftmost and rightmost bit

set(BUILTIN_POPCNT NO)
# Check if we are on a Linux system
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
	# Use /proc/cpuinfo to get the information
	file(STRINGS "/proc/cpuinfo" _cpuinfo)
	if(_cpuinfo MATCHES "(sse4_2)|(sse4a)")
		set(BUILTIN_POPCNT YES)
	endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
#  handle windows
#	get_filename_component(_vendor_id "[HKEY_LOCAL_MACHINE\\Hardware\\Description\\System\\CentralProcessor\\0;VendorIdentifier]" NAME CACHE)
#	get_filename_component(_cpu_id "[HKEY_LOCAL_MACHINE\\Hardware\\Description\\System\\CentralProcessor\\0;Identifier]" NAME CACHE)	
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
#  handle MacOs
endif()	
	
