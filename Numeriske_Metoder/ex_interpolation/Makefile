CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
CC = clang


.Phoney: default
default:	initial.data out.txt Makefile ex1_plot.svg ex2_qspline_const_plot.svg ex2_qspline_lin_plot.svg ex2_qspline_quad_plot.svg ex3_cspline_cos_plot.svg
	cat out.txt

out.txt plot.data:main initial.data
	./main  $$(wc --lines initial.data) < initial.data 1> out.txt  2> plot.data

ex1_plot.svg:plot.gpi plot.data initial.data
	gnuplot $< > $@

ex2_qspline_const_plot.svg:plot2.gpi qspline.data qinitial.data
	gnuplot $< > $@

ex2_qspline_lin_plot.svg:plot3.gpi qspline.data qinitial.data
	gnuplot $< > $@

ex2_qspline_quad_plot.svg:plot4.gpi qspline.data qinitial.data
		gnuplot $< > $@

ex3_cspline_cos_plot.svg:plot5.gpi cspline.data initial.data
		gnuplot $< > $@

qspline.o: qspline.c

cspline.o: cspline.c

myinteglin.o: myinteglin.c

main: main.o

initial.data: trigdata
	./$< > initial.data

trigdata: trigdata.o

clean:
	$(RM) initial.data out* *.o trigdata qspline.data ex3_cspline_cos_plot.svg cspline.data qinitial.data ex2_qspline_quad_plot.svg ex2_qspline_lin_plot.svg ex2_qspline_const_plot.svg ex1_plot.svg plot.data
