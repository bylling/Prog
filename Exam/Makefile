CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
CC = gcc

default: out.txt report.pdf Makefile
	cat $<

report.pdf: report.tex plot-cairo.tex plot-cairo-2.tex plot-cairo_time.tex
	pdflatex $<
	pdflatex $<

timeinverse_iter.txt:timeinverse_iter
	n=100; $(RM) $@; \
        while [ $${n} -le 800 ] ; do \
        	\time --format "%U" --append --output $@ ./$< $$n > /dev/null; \
        	n=`expr $$n + 100`; \
        done;

timejacobi.txt:timejacobi
	n=100; $(RM) $@; \
        while [ $${n} -le 800 ] ; do \
        	\time --format "%U" --append --output $@ ./$< $$n > /dev/null; \
        	n=`expr $$n + 100`; \
        done;

plot-cairo.tex: plot1.gpi plot1.data
	sed -e 's/TERM/cairolatex pdf/' -e 's/OUT/$@/' < $< | gnuplot

plot-cairo-2.tex: plot2.gpi plot2.data
	sed -e 's/TERM/cairolatex pdf/' -e 's/OUT/$@/' < $< | gnuplot

plot-cairo_time.tex: plot_time.gpi timeinverse_iter.txt timejacobi.txt
	sed -e 's/TERM/cairolatex pdf/' -e 's/OUT/$@/' < $< | gnuplot

plot2.data: figure2code
	./figure2code 2> error.txt

plot1.data: figure1code
	./figure1code 2> error.txt

out.txt: main
	./main 1> out.txt  2> error.txt

timeinverse_iter: timeinverse_iter.o functionfile.o function_lineq.o function_eigenval.o

timejacobi: timejacobi.o functionfile.o function_lineq.o function_eigenval.o

main:main.o functionfile.o function_lineq.o function_eigenval.o

figure1code:figure1code.o functionfile.o function_lineq.o function_eigenval.o

figure2code:figure2code.o functionfile.o function_lineq.o function_eigenval.o


clean:
	$(RM) *.o out* *.svg main *.pdf *.log *.aux *.data error.txt plot-cairo.tex plot-cairo_time.tex plot-cairo-2.tex timejacobi.txt timemain figure1code figure2code timejacobi timeinverse_iter timeinverse_iter.txt
