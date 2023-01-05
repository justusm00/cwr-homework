tests: my_numerics.o tests.o
	gcc -I.. -Wall my_numerics.o tests.o -lm -o  tests
my_numerics.o: ../my_numerics.c ../my_numerics.h
	gcc -c -I.. -Wall ../my_numerics.c
tests.o: my_numerics.o tests.c
	gcc -c -I.. -Wall tests.c
