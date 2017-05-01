hd= /Users/caojunzhi


CC= gcc
CLINK= gcc 
CFLAGS=  $(FFTW_CFLAGS)
CLIB= -L./ -lcutil -lm
XDIR= $(HOME)/exec


OBJS1=	main.o ran1.o gasdev.o zbrent.o erfinv.o get_redshift.o
merger_tree:	$(OBJS1)
	$(CLINK) $(CFLAGS) -o $@ $(OBJS1) $(CLIB)
	cp $@ $(HOME)/exec/

clean:
	rm -f *.o
