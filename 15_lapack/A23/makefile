a23: a23.o my_numerics.o
	gcc -Wall -Ofast -flto a23.o my_numerics.o -lgsl -lgslcblas -lm -o a23

a23.o: a23.c ../../my_numerics.h
	gcc -I../.. -Wall -Ofast -flto -c a23.c

my_numerics.o: ../../my_numerics.h ../../my_numerics.c
	gcc -Wall -Ofast -flto -c ../../my_numerics.c
