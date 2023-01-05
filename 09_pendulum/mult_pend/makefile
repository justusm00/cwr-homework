PENDULUMS: pendulums.o my_numerics.o
	gcc -Wall -Ofast -flto pendulums.o my_numerics.o -lgsl -lgslcblas -lm -o PENDULUMS

pendulums.o: pendulums.c ../../my_numerics.h
	gcc -I../.. -Wall -Ofast -flto -c pendulums.c

my_numerics.o: ../../my_numerics.h ../../my_numerics.c
	gcc -Wall -Ofast -flto -c ../../my_numerics.c
