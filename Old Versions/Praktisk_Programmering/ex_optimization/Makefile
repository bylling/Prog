CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)

default: out.txt Makefile plot.svg
	cat out.txt

out.txt plot.data:main Makefile rosen.o decayfunction.o
	./main $$(wc --lines decay.data) < decay.data 1> out.txt 2> plot.txt

plot.svg: plot.gpi plot.txt
	gnuplot $< >$@


clean:
	$(RM) main out* *.o plot.svg plot.txt
