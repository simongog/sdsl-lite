# Check if the CPU provides fast operations
# for popcount, leftmost and rightmost bit

set(BUILTIN_POPCNT 0)
if(DEFINED ENV{BUILD_PORTABLE})
	#message(STATUS "sse4.2 disabled")
# Check if we are on a Linux system
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
	# Use /proc/cpuinfo to get the information
	file(STRINGS "/proc/cpuinfo" _cpuinfo)
	if(_cpuinfo MATCHES "(sse4_2)|(sse4a)")
		set(BUILTIN_POPCNT 1)
	endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
#  handle windows
#	get_filename_component(_vendor_id "[HKEY_LOCAL_MACHINE\\Hardware\\Description\\System\\CentralProcessor\\0;VendorIdentifier]" NAME CACHE)
#	get_filename_component(_cpu_id "[HKEY_LOCAL_MACHINE\\Hardware\\Description\\System\\CentralProcessor\\0;Identifier]" NAME CACHE)	
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
#  handle MacOs
execute_process(COMMAND sysctl -n machdep.cpu.features
                OUTPUT_VARIABLE _cpuinfo OUTPUT_STRIP_TRAILING_WHITESPACE)
	if(_cpuinfo MATCHES "SSE4.2")
		set(BUILTIN_POPCNT 1)
	endif()
endif()	
	
