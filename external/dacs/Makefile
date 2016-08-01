CPP=gcc
//CXX=clang++

OBJECTS=obj/basics.o obj/bitrankw32int.o obj/directcodes.o obj/dacs.o
//CPPFLAGS=-Wall -g3 
#CPPFLAGS=-Wall -O9 -g -DNDEBUG 
CPPFLAGS=-Wall -O3 
DEST=.


all: $(OBJECTS)
	ar rvs libdacs.a $(OBJECTS)

obj/basics.o: src/basics.c
	$(CPP) $(CPPFLAGS) -c $< -o $@
obj/bitrankw32int.o: src/bitrankw32int.c
	$(CPP) $(CPPFLAGS) -c $< -o $@
obj/directcodes.o: src/directcodes.c
	$(CPP) $(CPPFLAGS) -c $< -o $@

obj/dacs.o: src/dacs.cc
	$(CXX) $(CPPFLAFS) -c $< -o $@


				
clean:
	rm -f $(OBJECTS)
	cd $(DEST); rm -f *.a

