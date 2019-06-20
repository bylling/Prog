CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
CC = gcc

default: out.txt report.pdf Makefile
	cat $<

report.pdf: report.tex plot-cairo.tex plot-cairo2.tex plot-cairo3.tex plot-cairo4.tex
	pdflatex $<
	pdflatex $<


plot-cairo.tex: plot.gpi plot.data
	sed -e 's/TERM/cairolatex pdf/' -e 's/OUT/$@/' < $< | gnuplot

plot-cairo2.tex: plot2.gpi plot.data
	sed -e 's/TERM/cairolatex pdf/' -e 's/OUT/$@/' < $< | gnuplot

plot-cairo3.tex: plot3.gpi plot.data
		sed -e 's/TERM/cairolatex pdf/' -e 's/OUT/$@/' < $< | gnuplot

plot-cairo4.tex: plot4.gpi plot.data
		sed -e 's/TERM/cairolatex pdf/' -e 's/OUT/$@/' < $< | gnuplot

out.txt plot.data:main
	./main -3 3 0.2  1> out.txt 2> plot.data

main:main.o myarctan.o myarccot.o myarctan2.o myarccot2.o

clean:
	$(RM) *.o out* *.svg main plot.data *.pdf *.log *.aux plot-cairo.tex plot-cairo2.tex plot-cairo3.tex plot-cairo4.tex
