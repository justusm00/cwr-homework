SINGLE_PEND: single_pend.o my_numerics.o
	gcc -Wall -Ofast -flto  single_pend.o my_numerics.o -lm -o SINGLE_PEND
	
	
single_pend.o: single_pend.c ../../my_numerics.h
	gcc -Wall -Ofast -flto  -I../.. -c single_pend.c

my_numerics.o: ../../my_numerics.h ../../my_numerics.c
	gcc -Wall -Ofast -flto -c ../../my_numerics.c 
