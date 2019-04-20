CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
CC = clang

.Phoney: default
default:	out.txt Makefile
	cat out.txt

out.txt: main
	./main 1> out.txt  2> error.txt

main: main.o rootfunctions.o lineqfunctions.o

clean:
	$(RM) out* *.o *.svg error.txt main
