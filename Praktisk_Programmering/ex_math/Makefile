CFLAGS = -std=gnu99 -Wall
LDLIBS = -lm
CC = gcc

.Phoney: default
default: out.txt
		cat out.txt

out.txt: main
	./main > out.txt

main: main.o
	$(CC) $(CFLAGS) $(LDFLAGS) $< $(LDLIBS) -o $@

main.o: main.c
	$(CC) $(CFLAGS) -c $<

.PHONEY: clean

clean:
	$(RM) main main.o out.txt
