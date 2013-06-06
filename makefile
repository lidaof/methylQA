include ../../inc/common.mk

L += $(MYSQLLIBS) -lm -lz
MYLIBDIR = ../../lib/${MACHTYPE}
MYLIBS =  ${MYLIBDIR}/jkhgap.a ${MYLIBDIR}/jkweb.a

MYF = methylQA

O = generic.o medip.o mre.o density.o genomecov.o from_kent.o $(MYF).o

all: ${O} $(MYLIBS)
	${CC} ${COPT} -o $(MYF) $O ${MYLIBS} $L
	cp $(MYF) ~/bin/
clean:
	rm -f methylQA *.o
