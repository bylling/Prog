CFLAGS = -Wall -std=c99 -O
LDLIBS = -lm
CC = gcc

.Phoney: default
default: out.txt
	 cat $<

out.txt: main
	./$< > $@

main: main.o nvector.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

main.o: main.c nvector.h
	$(CC) $(CFLAGS) -c $<

nvector.o: nvector.c nvector.h
	$(CC) $(CFLAGS) -c $<

.PHONEY: clean

clean:
	$(RM) main out* *.o *~

indent:
		indent --linux-style *.c; rm -f *.c~
