rootPath = ../
include ../include.mk

all: 
	cd lastz-distrib-1.03.54 && ${MAKE} install
	mv ${BINDIR}/lastz ${BINDIR}/cPecanLastz
	mv ${BINDIR}/lastz_D ${BINDIR}/cPecanLastz_D

clean: 
	cd lastz-distrib-1.03.54 && ${MAKE} clean 
	rm -rf ${BINDIR}/cPecanLastz ${BINDIR}/cPecanLastz_D

