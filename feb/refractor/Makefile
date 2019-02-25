# target: dependencies
#	action
#output: test.o simulator.o
#	g++-6 -fopenmp test.o simulator.o -o output

output: main.o simulator.o randomEngine.o
	g++-6 -fopenmp main.o simulator.o randomEngine.o -o output

#output: test_simulate.o simulator.o randomEngine.o
#	g++-6 -fopenmp test_simulate.o simulator.o randomEngine.o -o output


randomEngine.o : randomEngine.cpp RandomEngine.h Vector.h
	g++-6 -fopenmp -c randomEngine.cpp

test.o: test.cpp Vector.h simulator.h simulator.cpp
	g++-6 -fopenmp -c test.cpp

simulator.o: simulator.cpp simulator.h Vector.h randomEngine.h
	g++-6 -fopenmp -c simulator.cpp



main.o: main.cpp vector.h simulator.h simulator.cpp
	g++-6 -fopenmp -c main.cpp

test_simulate.o: test_simulate.cpp vector.h simulator.h simulator.cpp
	g++-6 -fopenmp -c test_simulate.cpp


clean:
	rm *.o
