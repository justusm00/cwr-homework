pend_self: my_numerics.o pend_self.o
	gcc -I.. -Wall -Ofast -flto  my_numerics.o pend_self.o -lm -o pend_self
my_numerics.o: ../my_numerics.c ../my_numerics.h
	gcc -c -I.. -Wall -Ofast -flto  ../my_numerics.c
pend_self.o: my_numerics.o pend_self.c
	gcc -c -I.. -Wall -Ofast -flto  pend_self.c
