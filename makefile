LINKLIB=  -L/home/samir/astro/lib  -fopenmp -lfftw3 -lcfitsio -lnrcp -lm -lgsl -lgslcblas
INCLUDE=-I/home/samir/astro/include/

corrl_v2_sim: corrl_v2_sim.c read_fits_func.c
	gcc -g -o corrl_v2_sim $(INCLUDE) read_fits_func.c corrl_v2_sim.c $(LINKLIB)
	rm -rf *~
