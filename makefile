KENT=/opt/kent/src
SAMTOOLS=/opt/samtools-0.1.17

CC=gcc
COPT= -O -g
CFLAGS= -Wall -Werror -Wformat -Wimplicit -Wreturn-type -Wuninitialized
DFLAGS= -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE
INCLUDES= -I$(KENT)/inc -I$(KENT)/hg/inc -I$(SAMTOOLS)
L += -lm -lz
MYLIBDIR = $(KENT)/lib/${MACHTYPE}
MYLIBS =  ${MYLIBDIR}/jkhgap.a ${MYLIBDIR}/jkweb.a

MYF = methylQA
O = generic.o medip.o mre.o density.o genomecov.o from_kent.o $(MYF).o

%.o: %.c
	${CC} ${COPT} ${CFLAGS} ${DFLAGS} $(INCLUDES) -o $@ -c $<


all: ${O} $(MYLIBS)
	${CC} ${COPT} -o $(MYF) $O ${MYLIBS} $(SAMTOOLS)/libbam.a -pthread -lssl $L
	cp $(MYF) ~/bin/
clean:
	rm -f methylQA *.o
