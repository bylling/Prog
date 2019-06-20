CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
CC = clang

.Phoney: default
default: out.txt Makefile plot.svg
		cat out.txt

plot.svg: plot.gpi plot.txt
		gnuplot $< > $@

main: main.o

out.txt plot.txt:main
		./$< 1> out.txt 2> plot.txt

clean:
		$(RM) main main.o out.txt plot.svg plot.txt
