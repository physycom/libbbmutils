test: obj/matrix_lib.o src/test.c
	gcc -o bin/rotations.exe obj/matrix_lib.o src/test.c
	./bin/rotations.exe

obj/matrix_lib.o: src/matrix_lib.h src/matrix_lib.c
	gcc -c -o obj/matrix_lib.o src/matrix_lib.c
