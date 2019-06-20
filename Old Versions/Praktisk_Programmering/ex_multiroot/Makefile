CFLAGS = -Wall -std=c11 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
CC = clang

blah: out.txt Makefile plot.svg
	cat out.txt

out.txt plot.txt: main
	./$< 1> out.txt   2> plot.txt

plot.svg: plot.gpi plot.txt
	gnuplot $< >$@


main: main.o root.o

clean:
	$(RM) main out* *.o plot.svg plot.txt
