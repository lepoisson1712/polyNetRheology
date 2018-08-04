CC      = /usr/bin/g++
CFLAGS  = -std=c++11 -Wall -pg -O3

#OBJ 	= polymer20160316.o defs.o simulation20160316.o crosslink.o crosslinkPair.o crosslinkPair2.o crosslinkTriple.o crosslinkTriple2_20160316.o rngmit.o
OBJ 	= polymer.o defs.o simulation.o crosslink.o  crosslinkPair.o  crosslinkTriple.o rngmit.o
#OBJ 	= polymer.o defs.o simulation_small.o crosslink.o  crosslinkPair2.o  crosslinkTriple2.o rngmit.o
#OBJ 	= polymer2.o defs.o simulation2.o crosslink.o crosslinkPair.o crosslinkPair2.o crosslinkTriple.o crosslinkTriple2.o rngmit.o


gransphero: $(OBJ)
	#$(CC) $(CFLAGS) -lboost_program_options -o polymer_oldPull $(OBJ)
	#$(CC) $(CFLAGS) -lboost_program_options -o polymer_newPull $(OBJ)
	#$(CC) $(CFLAGS) -lboost_program_options -o polymer_dummy $(OBJ)
	$(CC) $(CFLAGS) -lboost_program_options -o polymer $(OBJ)
	#$(CC) $(CFLAGS) -lboost_program_options -o polymer_lmin0.05smallTest $(OBJ)
%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

clean:
	rm *o
