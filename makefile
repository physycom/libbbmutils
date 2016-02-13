test: all
	./bin/test_3d.exe > bin/test_3d.log
	./bin/test_6d.exe > bin/test_6d.log

all: test_3d test_6d

test_3d: obj/math_lib.o src/test_3d.c
	gcc -o bin/test_3d.exe obj/math_lib.o src/test_3d.c

test_6d: obj/math_lib.o src/test_6d.c
	gcc -o bin/test_6d.exe obj/math_lib.o src/test_6d.c

obj/math_lib.o: src/math_lib.h src/math_lib.c
	gcc -c -o obj/math_lib.o src/math_lib.c

clean: 
	rm -rf bin/*.exe bin/*.log obj/math_lib.o