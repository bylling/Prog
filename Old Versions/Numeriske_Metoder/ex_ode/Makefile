CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
CC = gcc

.Phoney: default
default:	out.txt Makefile ex_A_plot.svg ex_B_plot.svg ex_C_plot.svg
	cat out.txt

ex_A_plot.svg:plotA.gpi plotdata.txt
	gnuplot $< > $@

ex_B_plot.svg:plotB.gpi plotdata.txt plotdata2.txt
		gnuplot $< > $@

ex_C_plot.svg:plotC.gpi plotdata3.txt
		gnuplot $< > $@

out.txt plotdata3.txt plotdata2.txt plotdata1.txt: main
	./main 1> out.txt  2> error.txt

main: main.o odefunctions.o integrationfunctions.o

clean:
	$(RM) out* *.o *.svg error.txt main plotdata.txt plotdata2.txt plotdata3.txt
