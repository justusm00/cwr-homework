LAPACK_LIBDIR =
LAPACK_INCDIR =
# You might need (and adapt) the following two lines if your lapacke
# libs and header are installed in a non-standard location
#LAPACK_LIBDIR = -L/usr/local/opt/lapack/lib
#LAPACK_INCDIR = -I/usr/local/opt/lapack/include

a22: a22.o
	gcc -Wall -Wextra -Wpedantic -Ofast -flto a22.o $(LAPACK_LIBDIR) -l lapacke -lm -o a22

a22.o: a22.c ../my_numerics.h
	gcc -Wall -Wextra -Wpedantic -I.. $(LAPACK_INCDIR) -c -flto -Ofast a22.c

run: a22
	./a22

clean:
	@rm -rf *.o
	@rm -f a22
