CFLAGS = -Wall -std=gnu1x $$(gsl-config --cflags)
LDLIBS = $$(gsl-config --libs)
CC = gcc


.Phoney: default
default: out.txt Makefile ex_A_plot.svg ex_B_plot.svg ex_C_plot.svg
	cat out.txt

out.txt plot.data:main
	./main  1> out.txt  2> error.txt

ex_A_plot.svg: plot.gpi time_cyclic.txt
	gnuplot $< > $@

ex_B_plot.svg: plot_2.gpi time_cyclic.txt time_mod_1.txt time_mod_all.txt
		gnuplot $< > $@

ex_C_plot.svg: plot_3.gpi time_cyclic.txt time_c.txt
		gnuplot $< > $@

time_cyclic.txt: jacobi_cyclic_time
	n=100; $(RM) $@; \
	        while [ $${n} -le 500 ] ; do \
	        	\time --format "%U" --append --output $@ ./$< $$n > /dev/null; \
	        	n=`expr $$n + 30`; \
	        done;

time_mod_1.txt: jacobi_mod_time
	n=100; $(RM) $@; \
	        while [ $${n} -le 500 ] ; do \
	        	\time --format "%U" --append --output $@ ./$< $$n > /dev/null; \
	        	n=`expr $$n + 30`; \
	        done;

time_mod_all.txt: jacobi_mod_all_time
	n=100; $(RM) $@; \
	        while [ $${n} -le 500 ] ; do \
	        	\time --format "%U" --append --output $@ ./$< $$n > /dev/null; \
	        	n=`expr $$n + 30`; \
	        done;

time_c.txt: jacobi_time_c
		n=100; $(RM) $@; \
	        while [ $${n} -le 300 ] ; do \
	        	\time --format "%U" --append --output $@ ./$< $$n > /dev/null; \
	        	n=`expr $$n + 30`; \
	        done;

jacobi_cyclic_time: jacobi_cyclic_time.o eigenfunctions.o

jacobi_mod_all_time: jacobi_mod_all_time.o eigenfunctions.o

jacobi_mod_time: jacobi_mod_time.o eigenfunctions.o

jacobi_time_c: jacobi_time_c.o eigenfunctions.o

main: main.o eigenfunctions.o

clean:
	$(RM) out* *.o main jacobi_cyclic_time time_cyclic.txt ex_A_plot.svg time_mod_1.txt ex_B_plot.svg jacobi_mod_time time_mod_all.txt jacobi_mod_all_time time_c.txt jacobi_time_c ex_C_plot.svg error.txt
