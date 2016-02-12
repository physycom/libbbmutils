test: obj/matrix_lib.o src/test.c
	gcc -o bin/test.exe obj/matrix_lib.o src/test.c
	./bin/test.exe

obj/matrix_lib.o: src/matrix_lib.h src/matrix_lib.c
	gcc -c -o obj/matrix_lib.o src/matrix_lib.c

clean: 
	rm -rf bin/test.exe obj/matrix_lib.o