CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
CC = gcc


.Phoney: default
default:	initial.data out.txt Makefile ex_A_plot.svg ex_B_qspline_const_plot.svg ex_B_qspline_lin_plot.svg ex_B_qspline_quad_plot.svg ex_C_cspline_cos_plot.svg
	cat out.txt

out.txt plot.data:main initial.data
	./main  $$(wc --lines initial.data) < initial.data 1> out.txt  2> error.txt

ex_A_plot.svg:plot.gpi plot.data initial.data
	gnuplot $< > $@

ex_B_qspline_const_plot.svg:plot2.gpi qspline.data qinitial.data
	gnuplot $< > $@

ex_B_qspline_lin_plot.svg:plot3.gpi qspline.data qinitial.data
	gnuplot $< > $@

ex_B_qspline_quad_plot.svg:plot4.gpi qspline.data qinitial.data
		gnuplot $< > $@

ex_C_cspline_cos_plot.svg:plot5.gpi cspline.data initial.data
		gnuplot $< > $@


main: main.o interpfunctions.o

initial.data: trigdata
	./$< > initial.data

trigdata: trigdata.o

clean:
	$(RM) initial.data out* *.o *.svg trigdata qspline.data cspline.data qinitial.data plot.data
