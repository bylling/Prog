CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = `gsl-config --libs`

all: airy.svg out.txt
	cat out.txt

airy.svg: airy.gpi airy.txt
	gnuplot $< > $@

out.txt airy.txt:main
	cat input.matrix | ./$< 1> out.txt 2> airy.txt

main: main.o make_airy_plot.o

main.o make_airy_plot.o: make_airy_plot.h

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)

.PHONEY:clean
clean:
	$(RM) *.o out* *.svg main airy.txt
