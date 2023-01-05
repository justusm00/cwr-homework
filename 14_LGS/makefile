test: test.o my_numerics.o
	gcc -Wall test.o my_numerics.o -lm -o test

test.o: test.c ../my_numerics.h
	gcc -I.. -Wall -c test.c

my_numerics.o: ../my_numerics.h ../my_numerics.c
	gcc -Wall -c ../my_numerics.c
