max: max.o my_numerics.o
	gcc -Wall -Ofast -flto max.o my_numerics.o -lgsl -lgslcblas -lm -o max

max.o: max.c ../../my_numerics.h
	gcc -I../.. -Wall -Ofast -flto -c max.c

my_numerics.o: ../../my_numerics.h ../../my_numerics.c
	gcc -Wall -Ofast -flto -c ../../my_numerics.c
