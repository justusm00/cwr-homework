derivative: my_numerics.o derivative.o
	gcc -I.. -Wall my_numerics.o derivative.o -lm -o  derivative
my_numerics.o: ../my_numerics.c ../my_numerics.h
	gcc -c -I.. -Wall ../my_numerics.c
derivative.o: my_numerics.o derivative.c
	gcc -c -I.. -Wall derivative.c
