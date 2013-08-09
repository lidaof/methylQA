CC=gcc
COPT= -O -g
CFLAGS= -Wall -Werror -Wformat -Wimplicit -Wreturn-type -Wuninitialized
DFLAGS= -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE
KENT=cuskent
SAMTOOLS=cussamtools
INCLUDES= -I$(KENT) -I$(SAMTOOLS)
KENTLIB=$(KENT)/libcuskent.a
SAMLIB=$(SAMTOOLS)/libbam.a
L = -pthread -lm -lz

MYF = methylQA
O = generic.o medip.o mre.o density.o genomecov.o from_kent.o $(MYF).o

%.o: %.c
	$(CC) $(COPT) $(CFLAGS) $(DFLAGS) $(INCLUDES) -o $@ -c $<

all: $O $(KENTLIB) $(SAMLIB)
	$(CC) $(COPT) -o $(MYF) $O $(KENTLIB) $(SAMLIB) $L

.PHONY:all $(KENTLIB) $(SAMLIB) clean

$(KENTLIB):
	cd $(KENT) && make

$(SAMLIB):
	cd $(SAMTOOLS) && make lib

cleanlocal:
	rm -f $(MYF) $(O)

clean:
	wdir=`pwd`; \
	cd $$wdir/$(KENT) && make clean; \
	cd $$wdir/$(SAMTOOLS) && make cleanlocal; \
	cd $$wdir && rm -f $(MYF) $(O)
