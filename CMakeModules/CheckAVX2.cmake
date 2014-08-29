# Check if the CPU provides fast operations
# for popcount, leftmost and rightmost bit

set(AVX2 0)
# Check if we are on a Linux system
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
	# Use /proc/cpuinfo to get the information
	file(STRINGS "/proc/cpuinfo" _cpuinfo)
	if(_cpuinfo MATCHES "(avx2)")
    set(AVX2 1)
	endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
#  handle windows
#	get_filename_component(_vendor_id "[HKEY_LOCAL_MACHINE\\Hardware\\Description\\System\\CentralProcessor\\0;VendorIdentifier]" NAME CACHE)
#	get_filename_component(_cpu_id "[HKEY_LOCAL_MACHINE\\Hardware\\Description\\System\\CentralProcessor\\0;Identifier]" NAME CACHE)	
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
#  handle MacOs
execute_process(COMMAND sysctl -n machdep.cpu.features
                OUTPUT_VARIABLE _cpuinfo OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(_cpuinfo MATCHES "AVX2")
    set(AVX2 1)
	endif()
endif()	
	
