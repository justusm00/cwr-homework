lagrange: my_numerics.o lagrange.o
	gcc -I.. -Wall my_numerics.o lagrange.o -lm -o  lagrange
my_numerics.o: ../my_numerics.c ../my_numerics.h
	gcc -c -I.. -Wall ../my_numerics.c
lagrange.o: my_numerics.o lagrange.c
	gcc -c -I.. -Wall lagrange.c
