CPP=gcc

OBJECTS=basic.o bitrankw32int.o directcodes.o kTree.o hash.o MemoryManager.o
BINS=build_tree  test_tree rebuild_tree invrebuild_tree revtest_tree compress_leaves checkrebuild_tree rangerebuild_tree use_tree speed_direct speed_reverse

#CPPFLAGS=-Wall -g3 
CPPFLAGS=-O9 -g -DNDEBUG 
DEST=.

%.o: %.c
	$(CPP) $(CPPFLAGS) -c $< -o $@

all: bin

bin: $(OBJECTS) $(BINS)

build_tree:
	gcc $(CPPFLAGS) -o $(DEST)/build_tree buildtree.c $(OBJECTS) -lm

compress_leaves:
	gcc $(CPPFLAGS) -o $(DEST)/compress_leaves compressleaves.c $(OBJECTS) -lm

use_tree:
	gcc $(CPPFLAGS) -o $(DEST)/use_tree usetree.c $(OBJECTS) -lm

test_tree:
	gcc $(CPPFLAGS) -o $(DEST)/test_tree speedtree.c $(OBJECTS) -lm

revtest_tree:
	gcc $(CPPFLAGS) -o $(DEST)/revtest_tree speedreverse.c $(OBJECTS) -lm
	
rebuild_tree:
	gcc $(CPPFLAGS) -o $(DEST)/rebuild_tree rebuildtree.c $(OBJECTS) -lm

invrebuild_tree:
	gcc $(CPPFLAGS) -o $(DEST)/invrebuild_tree rebuildreverse.c $(OBJECTS) -lm

checkrebuild_tree:
	gcc $(CPPFLAGS) -o $(DEST)/checkrebuild_tree rebuildchecklink.c $(OBJECTS) -lm

rangerebuild_tree:
	gcc $(CPPFLAGS) -o $(DEST)/rangerebuild_tree rebuildrange.c $(OBJECTS) -lm

speed_direct:
	gcc $(CPPFLAGS) -o $(DEST)/speed_direct speedtree.c $(OBJECTS) -lm

speed_reverse:
	gcc $(CPPFLAGS) -o $(DEST)/speed_reverse speedreverse.c $(OBJECTS) -lm

				
clean:
	rm -f $(OBJECTS) $(BINS)
	cd $(DEST); rm -f *.a $(BINS)

