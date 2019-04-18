CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
CC = clang

.Phoney: default
default:	out.txt Makefile ex_A_plot.svg ex_B_plot.svg ex_C_plot.svg
	cat out.txt

out.txt data1.txt data2.txt:main initial.data
	./main  $$(wc --lines initial.data) < initial.data 1> out.txt  2> error.txt

ex_A_plot.svg:plot.gpi initial.data data1.txt
	gnuplot $< > $@

ex_B_plot.svg:plot2.gpi initial.data data2.txt
	gnuplot $< > $@

ex_C_plot.svg:plot3.gpi initial.data data3.txt
	gnuplot $< > $@


main: main.o least_square_functions.o lineqfunctions.o eigenfunctions.o

clean:
	$(RM) out* *.o *.svg error.txt data1.txt main data2.txt data3.txt
