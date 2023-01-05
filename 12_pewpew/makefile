COULOMB: coulomb.o my_numerics.o
	gcc -Wall coulomb.o my_numerics.o -lm -o COULOMB

coulomb.o: coulomb.c ../my_numerics.h
	gcc -I.. -Wall -c coulomb.c

my_numerics.o: ../my_numerics.h ../my_numerics.c
	gcc -Wall -c ../my_numerics.c
