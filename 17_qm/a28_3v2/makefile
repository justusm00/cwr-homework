LAPACK_LIBDIR =
LAPACK_INCDIR =
# You might need (and adapt) the following two lines if your lapacke
# libs and header are installed in a non-standard location
#LAPACK_LIBDIR = -L/usr/local/opt/lapack/lib
#LAPACK_INCDIR = -I/usr/local/opt/lapack/include

a28_3: a28_3.o
	gcc -Wall -Wextra -Wpedantic -Ofast -flto a28_3.o $(LAPACK_LIBDIR) -l lapacke -lm -o a28_3

a28_3.o: a28_3.c ../../my_numerics.h
	gcc -Wall -Wextra -Wpedantic -I.. $(LAPACK_INCDIR) -c -flto -Ofast a28_3.c

run: a28_3
	./a28_3

clean:
	@rm -rf *.o
	@rm -f a28_3
