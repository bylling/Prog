CFLAGS = -pthread -std=gnu1x
LDLIBS = -pthread
CC = clang
all:out.txt
	cat $<

out.txt:main
	time --append --portability -o $@ ./$< > $@
#	time --append --format '%U %S %E' -o $@ ./$< > $@

clean:
	$(RM) main out.txt
