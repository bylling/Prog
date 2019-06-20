CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
CC = gcc

.Phoney: default
default:	out.txt Makefile ex_B_plot.svg
	cat out.txt

out.txt plotdata1.txt: main
	./main 1> out.txt  2> error.txt


ex_B_plot.svg: plot.gpi plotdata1.txt
	gnuplot $< > $@

main: main.o montecarlofunctions.o

clean:
	$(RM) out* *.o *.svg error.txt main plotdata1.txt
