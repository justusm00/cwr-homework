cannon: cannon.o my_numerics.o
	gcc -Wall -Ofast -flto cannon.o my_numerics.o -lgsl -lgslcblas -lm -o cannon

cannon.o: cannon.c ../../my_numerics.h
	gcc -I../.. -Wall -Ofast -flto -c cannon.c

my_numerics.o: ../../my_numerics.h ../../my_numerics.c
	gcc -Wall -Ofast -flto -c ../../my_numerics.c
