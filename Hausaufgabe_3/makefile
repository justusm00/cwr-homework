toeplitz: toeplitz.o my_numerics.o
	gcc -Wall toeplitz.o my_numerics.o -lm -o toeplitz

toeplitz.o: toeplitz.c ../my_numerics.h
	gcc -I.. -Wall -c toeplitz.c

my_numerics.o: ../my_numerics.h ../my_numerics.c
	gcc -Wall -c ../my_numerics.c
