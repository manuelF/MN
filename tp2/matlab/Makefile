BIN=../src/DCT


all: runfiles
	gnuplot *.plt  >&2 2> /dev/null

extra: latest
	$(BIN) 1 < ../datos/brian_kernighan.pgm

runfiles: latest 
	touch mod recovered orig
	$(BIN) < ../datos/dopp512.txt 
	mv -f mod mod1
	mv -f recovered rec1
	mv -f orig orig1
	touch mod recovered orig
	$(BIN) < ../datos/ramp1234.txt
	mv -f mod mod2
	mv -f recovered rec2
	mv -f orig orig2
	touch mod recovered orig
	$(BIN) < ../datos/dopp1024.txt
	mv -f mod mod3
	mv -f recovered rec3
	mv -f orig orig3
	touch mod recovered orig
	$(BIN) < ../datos/g450.txt
	mv -f mod mod4
	mv -f recovered rec4
	mv -f orig orig4
	touch mod recovered orig
	$(BIN) < ../datos/s64.txt
	mv -f mod mod5
	mv -f recovered rec5
	mv -f orig orig5


latest:
	cd ../src; make -f ../src/Makefile 

clean:
	rm mod* rec* orig*

