test: obj/math_lib.o src/test.c
	gcc -o bin/test.exe obj/math_lib.o src/test.c
	./bin/test.exe

obj/math_lib.o: src/math_lib.h src/math_lib.c
	gcc -c -o obj/math_lib.o src/math_lib.c

clean: 
	rm -rf bin/test.exe obj/math_lib.o