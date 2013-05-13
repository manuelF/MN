BIN=../src/DCT

all: runfiles
	gnuplot *.plt

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

latest:
	cd ../src; make -f ../src/Makefile
