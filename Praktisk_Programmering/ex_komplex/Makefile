CFLAGS = -Wall -std=gnu99 -O
LDLIBS = -lm
CC = gcc

.Phoney: default
default: out.txt
	 cat $<

out.txt: main
	./$< > $@

main: main.o komplex.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

main.o: main.c komplex.h
	$(CC) $(CFLAGS) -c $<

komplex.o: komplex.c komplex.h
	$(CC) $(CFLAGS) -c $<

.PHONEY: clean

clean:
	$(RM) main out* *.o *~

STYLE = --linux-style
.PHONEY:indent
indent           :                  ; indent $(STYLE) main.c komplex.c; rm -f *.c~
