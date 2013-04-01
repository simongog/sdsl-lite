LIB_DIR = @CMAKE_INSTALL_PREFIX@/lib
INC_DIR = @CMAKE_INSTALL_PREFIX@/include

# Returns $1-th .-separated part of string $2.
dim = $(word $1, $(subst ., ,$2))

# Returns value stored in column $3 for item with ID $2 in 
# config file $1
config_select=$(shell cat $1 | grep -v "^\#" | grep "$2;" | cut -f $3 -d';' )

# Returns value stored in column $3 for a line matching $2
# in config file $1
config_filter=$(shell cat $1 | grep -v "^\#" | fgrep "$2" | cut -f $3 -d';' )

# Get all IDs from a config file $1
config_ids=$(shell cat $1 | grep -v "^\#" | cut -f 1 -d';')

# Get column $2 from a config file $1
config_column=$(shell cat $1 | grep -v "^\#" | cut -f $2 -d';')

# Get size of file $1 in bytes
file_size=$(shell wc -c < $1 | tr -d ' ')

../data/%: 
	$(eval URL:=$(call config_filter,test_case.config,$@,4))
	@$(if $(URL),,\
		$(error "No downlaod link nor generation program specified for test case $@") )
	@echo "Downlaod input from $(URL) using curl"
	$(eval DEST_DIR:=$(shell dirname $@))
	cd $(DEST_DIR); curl -O $(URL)
	$(eval FILE:=$(DEST_DIR)/$(notdir $(URL)))
	@$(if $(filter-out ".gz",$(FILE)),\
		echo "Extract file $(FILE) using gunzip";\
		gunzip $(FILE))



