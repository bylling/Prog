CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
CC = clang

default: out.txt report.pdf Makefile
	cat $<

report.pdf: report.tex plot-cairo.tex
	pdflatex $<
	pdflatex $<


plot-cairo.tex: plot.gpi plot.data
	sed -e 's/TERM/cairolatex pdf/' -e 's/OUT/$@/' < $< | gnuplot

out.txt plot.data:main
	./main -3 3 0.2  1> out.txt 2> plot.data


main:main.o myerr.o

clean:
	$(RM) *.o out* *.svg main plot.data *.pdf *.log *.aux plot-cairo.tex
