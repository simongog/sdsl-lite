LIB_DIR = @CMAKE_INSTALL_PREFIX@/lib
INC_DIR = @CMAKE_INSTALL_PREFIX@/include

# Returns $1-th .-separated part of string $2.
dim = $(word $1, $(subst ., ,$2))

# Returns value stored in column $3 for item with ID $2 in 
# config file $1
config_select=$(shell cat $1 | grep "$2;" | cut -f $3 -d';' )

# Get all IDs from a config file $1
config_ids=$(shell cat $1 | grep -v "^\#" | cut -f 1 -d';')

