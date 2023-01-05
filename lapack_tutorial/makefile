LAPACK_LIBDIR =
LAPACK_INCDIR =
# You might need (and adapt) the following two lines if your lapacke
# libs and header are installed in a non-standard location
#LAPACK_LIBDIR = -L/usr/local/opt/lapack/lib
#LAPACK_INCDIR = -I/usr/local/opt/lapack/include

TUTORIAL: lapack_tutorial.o
	gcc -Wall -Wextra -Wpedantic -Ofast -flto lapack_tutorial.o $(LAPACK_LIBDIR) -l lapacke -lm -o TUTORIAL

lapack_tutorial.o: lapack_tutorial.c ../my_numerics.h
	gcc -Wall -Wextra -Wpedantic -I.. $(LAPACK_INCDIR) -c -flto -Ofast lapack_tutorial.c

run: TUTORIAL
	./TUTORIAL

clean:
	@rm -rf *.o
	@rm -f TUTORIAL
