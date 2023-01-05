POPULATION: population.o my_numerics.o
	gcc -Wall population.o my_numerics.o -lm -o POPULATION

population.o: population.c ../my_numerics.h
	gcc -I.. -Wall -c population.c

my_numerics.o: ../my_numerics.h ../my_numerics.c
	gcc -Wall -c ../my_numerics.c
