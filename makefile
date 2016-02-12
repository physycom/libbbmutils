all: matrix_lib.o rotations.c
	gcc -o rotations.exe matrix_lib.o rotations.c
	./rotations.exe

matrix_lib.o: matrix_lib.h matrix_lib.c
	gcc -c -o matrix_lib.o matrix_lib.c
