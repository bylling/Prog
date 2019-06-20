CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
CC = gcc

.Phoney: default
default:	out.txt Makefile ex_B_plot.svg
	cat out.txt

ex_B_plot.svg:plot.gpi data1.txt fitdata.txt
		gnuplot $< > $@

out.txt data1.txt fitdata.txt: main
	./main 1> out.txt  2> error.txt

main: main.o minimizationfunctions.o lineqfunctions.o rootfunctions.o downhillsimplex.o

clean:
	$(RM) out* *.o *.svg error.txt main ex_B_plot.svg data1.txt fitdata.txt
