CFLAGS = -Wall -std=gnu99 $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
CC = clang

default:out.txt plot.svg
	cat $<

plot.svg:plot.gpi plot.data
	gnuplot $< > $@

out.txt plot.data:main
	./main 1> out.txt 2> plot.data


main:main.o myorbit.o

test:
	echo $(PATH)

clean:
	$(RM) *.o out* *.svg main plot.data

INDENT_FLAGS = --linux-style
indent:
	for f in *.c; do indent $(INDENT_FLAGS) $$f; done
	rm -f *~
