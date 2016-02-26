measurements.o: 	measurements.c measurements.h
	gcc -O3 -c -fpic -o $@ measurements.c

measurements_wrap.o:	measurements_wrap.c
	gcc -O3 -c -fpic -o $@ measurements_wrap.c

measurements_wrap.c:	measurements.c measurements.h measurements.i
	swig -tcl8 measurements.i

cfacv.so:  cfacv.o cfacv_wrap.o measurements.o centers.o 
	gcc -O3 -shared -o $@ cfacv.o cfacv_wrap.o measurements.o centers.o -lm
	cp $@ ../lib

libgenericdataspace.so: genericdataspace.o genericdataspace_wrap.o
	gcc -O3 -shared -o $@ genericdataspace.o genericdataspace_wrap.o -lm
	cp $@ ../lib

genericdataspace.o:     genericdataspace.c genericdataspace.h
	gcc -O3 -c -fpic -o $@ genericdataspace.c

genericdataspace_wrap.o:        genericdataspace_wrap.c
	gcc -O3 -c -fpic -o $@ $?

genericdataspace_wrap.c:        genericdataspace.c genericdataspace.h genericdataspace.i
	swig -tcl8 genericdataspace.i

cfacv.o:      cfacv.c cfacv.h wrapcoords.h
	gcc -O3 -c -fpic -o $@ cfacv.c

cfacv_wrap.o: cfacv_wrap.c
	gcc -O3 -c -fpic -o $@ $?

cfacv_wrap.c: cfacv.c cfacv.i cfacv.h
	swig -tcl8 cfacv.i

centers.o:      centers.c centers.h
	gcc -O3 -c -fpic -o $@ centers.c

clean:
	rm *.o ../lib/*.so *~
