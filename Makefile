all: main

foralloperations.o: foralloperations.cpp foralloperations.hpp
	g++ -std=c++14 -Wall -g -c foralloperations.cpp

mesh.o: mesh.cpp mesh.hpp
	g++ -std=c++14 -Wall -g -c mesh.cpp

solution.o: solution.cpp solution.hpp
	g++ -std=c++14 -Wall -g -c solution.cpp

fields.o: fields.cpp fields.hpp mesh.hpp solution.hpp foralloperations.hpp
	g++ -std=c++14 -Wall -g -c fields.cpp

finitematrix.o: finitematrix.cpp finitematrix.hpp fields.hpp foralloperations.hpp
	g++ -std=c++14 -Wall -g -c finitematrix.cpp

equation.o: equation.cpp equation.hpp mesh.hpp fields.hpp finitematrix.hpp
	g++ -std=c++14 -Wall -g -c equation.cpp

filewrite.o: filewrite.cpp filewrite.hpp
	g++ -std=c++14 -Wall -g -c filewrite.cpp

main.o: main.cpp fields.hpp foralloperations.hpp initializevariables.hpp initializefinitematrixvar.hpp finitevolumeoperations.hpp filewrite.hpp
	g++ -std=c++14 -Wall -g -c main.cpp

main: main.o fields.o foralloperations.o solution.o mesh.o finitematrix.o equation.o filewrite.o
	g++ -floop-parallelize-all -ftree-parallelize-loops=8 -std=c++14 -Wall -g -lm -o main main.o fields.o foralloperations.o solution.o mesh.o finitematrix.o equation.o filewrite.o

clean:
	rm main *.o *.dat

