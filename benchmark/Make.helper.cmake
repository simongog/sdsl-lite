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

../data/%.z.info:
	$(eval TC:=../data/$*)
	@echo "Get xz-compression ratio for $(TC)"
	$(eval TC_XZ:=$(TC).xz)
	$(shell xz -9 -z -k -c $(TC) > $(TC_XZ))
	$(eval XZ_SIZE:=$(call file_size,$(TC_XZ)))
	$(shell rm $(TC_XZ)) 
	@echo "Get gzip-compression ratio for $(TC)"
	$(eval TC_GZ:=$(TC).gz) 
	$(shell gzip -9 -c $(TC) > $(TC_GZ))
	$(eval GZ_SIZE:=$(call file_size,$(TC_GZ)))
	$(shell rm $(TC_GZ))
	$(eval SIZE:=$(call file_size,$(TC)))
	$(eval XZ_RATIO:=$(shell echo "scale=2;100*$(XZ_SIZE)/$(SIZE)" | bc -q))
	$(eval GZ_RATIO:=$(shell echo "scale=2;100*$(GZ_SIZE)/$(SIZE)" | bc -q))
	@echo "xz;$(XZ_RATIO);xz -9\ngzip;$(GZ_RATIO);gzip -9" > $@
	


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


