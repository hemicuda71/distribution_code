all: rdmpi

rdmpi: read_dist_mpi.o alglibinternal.o ap.o specialfunctions.o timer.o 
	mpic++ read_dist_mpi.o alglibinternal.o ap.o specialfunctions.o timer.o -o rdmpi

read_dist_mpi.o:
	mpic++ -c -O3 read_dist_mpi.cpp

alglibinternal.o:
	g++ -c -O3 ./special_funcs/alglibinternal.cpp

ap.o:
	g++ -c -O3 ./special_funcs/ap.cpp

specialfunctions.o:
	g++ -c -O3 ./special_funcs/specialfunctions.cpp

timer.o:
	g++ -c -O3 timer.cpp

clean:
	rm rdmpi
	rm read_dist_mpi.o


