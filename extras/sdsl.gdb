#                                                                                                        
#   SDSL data structure inspector for GDB 
#
#  Inspired by the STL equivalent from
#  http://www.yolinux.com/TUTORIALS/src/dbinit_stl_views-1.03.txt
#
#   Copy the content of this file in your .gdbinit file to
#   use the functions.
#
#   The following SDSL containers are currently supported:
#
#       sdsl::int_vector<w> -- via the pv command

#
# sdsl::int_vector
#

define pv
	if $argc == 0
		help pv
	else
		set $width= $arg0.m_width
		set $size = $arg0.m_size / $width
	end
	set $i = 0
	set $end = $size-1
	if $argc == 2
		set $i = $arg1
		set $end = $arg1
	end
	if $argc == 3
		set $i = $arg1
		set $end = $arg2
		if $end > $size-1
			printf "ERROR: end=%u > %u=size-1. Please specify a valid range.", $end, $size-1
			set $i=1
			set $end=1
		end
	end	
	while $i <= $end
		set $offset = ($i*$width)%64
		set $word   = ($arg0.m_data+(($i*$width)/64))
		set $value = 0
		if $width > 0
			set $w1=(*$word)>>$offset
			if $offset+$width > 64
				set $s=64-(($offset+$width)&0x3F)
				set $value = $w1 | ((((*($word+1))<<$s)>>$s)<<(64-$offset))
			else
				set $s=64-$width  	
				set $value = ($w1<<$s)>>$s
			end
		end
		printf "elem[%u]: %u\n ", $i, $value
		set $i++
	end
	if $argc > 0 
		printf "int_vector size = %u\n", $size
	end	
end

document pv
	Prints std::int_vector<width> information.
	Syntax: pv <int_vector> <idx1> <idx2>
	Note: idx, idx1 and idx2 must be in the range [0..<int_vector>.size()-1].
	Examples:
	pv v - Prints vector content, size, capacity and T typedef
	pv v 0 - Prints element[idx] from int_vector
	pv v 1 2 - Prints elements in range [idx1..idx2] from int_vector
end
