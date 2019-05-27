CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
CC = gcc

.Phoney: default
default:	out.txt Makefile ex_A_plot.svg ex_B_plot.svg ex_C_plot.svg
	cat out.txt

out.txt plotdata1.txt plotdata2.txt: main
	./main 1> out.txt  2> error.txt

ex_A_plot.svg: plot.gpi plotdata1.txt
	gnuplot $< > $@

ex_B_plot.svg: plot2.gpi plotdata2.txt
	gnuplot $< > $@

ex_C_plot.svg: plot3.gpi plotdata3.txt
	gnuplot $< > $@

main: main.o neuralnetworkfunctions.o

clean:
	$(RM) out* *.o *.svg error.txt main plotdata1.txt plotdata2.txt plotdata3.txt
