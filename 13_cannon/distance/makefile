distance: distance.o my_numerics.o
	gcc -Wall -Ofast -flto distance.o my_numerics.o -lgsl -lgslcblas -lm -o distance

distance.o: distance.c ../../my_numerics.h
	gcc -I../.. -Wall -Ofast -flto -c distance.c

my_numerics.o: ../../my_numerics.h ../../my_numerics.c
	gcc -Wall -Ofast -flto -c ../../my_numerics.c
