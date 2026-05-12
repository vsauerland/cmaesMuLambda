CCC = g++
CC11 = gnu++0x
EIGEN = -I eigen-3.4.0/Eigen/
CCFLAGS = -O -std=$(CC11) 

PROGRAMS = testfunctions cmaes 

all:
	make $(PROGRAMS)

clean :
	/bin/rm -rf $(PROGRAMS)
	/bin/rm -rf *.out *.o *~ *.txt *.log

cmaes: auxiliaries.o cmaes.o
	$(CCC) $(CCFLAGS) $(EIGEN) auxiliaries.o cmaes.cpp -o cmaes
cmaes.o: cmaes.cpp
	$(CCC) -c $(CCFLAGS) $(EIGEN) cmaes.cpp -o cmaes.o
testfunctions: testfunctions.o
	$(CCC) $(CCFLAGS) testfunctions.o -o testfunctions
testfunctions.o: testfunctions.cpp
	$(CCC) -c $(CCFLAGS) $(EIGEN) testfunctions.cpp -o testfunctions.o
auxiliaries.o: auxiliaries.cpp
	$(CCC) -c $(CCFLAGS) $(EIGEN) auxiliaries.cpp -o auxiliaries.o
