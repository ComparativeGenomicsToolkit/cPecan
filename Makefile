rootPath = ./
include ./include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

cPecanDependencies =  ${LIBDEPENDS}
cPecanLibs = ${LDLIBS}
CFLAGS += -Wno-overlength-strings

all : ${LIBDIR}/cPecanLib.a ${BINDIR}/cPecanLibTests ${BINDIR}/cPecanRealign ${BINDIR}/cPecanEm ${BINDIR}/cPecanModifyHmm ${BINDIR}/cPecanAlign
	cd externalTools && ${MAKE} all

clean : 
	rm -f ${BINDIR}/cPecanRealign ${BINDIR}/cPecanEm ${BINDIR}/cPecanLibTests  ${LIBDIR}/cPecanLib.a
	cd externalTools && ${MAKE} clean
	rm -rf tmp_* *.o

export PYTHONPATH = ../sonLib/src:..
export PATH := ../sonLib/bin:${PATH}
test : all
	${PYTHON} allTests.py

${BINDIR}/cPecanRealign : cPecanRealign.c ${LIBDIR}/cPecanLib.a ${cPecanDependencies} 
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cPecanRealign cPecanRealign.c ${LIBDIR}/cPecanLib.a ${cPecanLibs} ${LDLIBS}

${BINDIR}/cPecanEm : cPecanEm.py
	cp cPecanEm.py ${BINDIR}/cPecanEm
	chmod +x ${BINDIR}/cPecanEm

${BINDIR}/cPecanModifyHmm : cPecanModifyHmm.py
	cp cPecanModifyHmm.py ${BINDIR}/cPecanModifyHmm
	chmod +x ${BINDIR}/cPecanModifyHmm

${BINDIR}/cPecanAlign : cPecanAlign.c ${LIBDIR}/cPecanLib.a ${cPecanDependencies} 
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -o ${BINDIR}/cPecanAlign cPecanAlign.c ${LIBDIR}/cPecanLib.a ${cPecanLibs} ${LDLIBS}

${BINDIR}/cPecanLibTests : ${libTests} ${LIBDIR}/cPecanLib.a ${cPecanDependencies}
	${CC} ${CPPFLAGS} ${CFLAGS} ${LDFLAGS} -Wno-error -o ${BINDIR}/cPecanLibTests ${libTests} ${LIBDIR}/cPecanLib.a ${cPecanLibs} ${LDLIBS}

${LIBDIR}/cPecanLib.a : ${libSources} ${libHeaders} ${stBarDependencies}
	${CC} ${CPPFLAGS} ${CFLAGS} -I inc -I ${LIBDIR}/ -c ${libSources} 
	${AR} rc cPecanLib.a *.o
	${RANLIB} cPecanLib.a 
	mv cPecanLib.a ${LIBDIR}/
	cp ${libHeaders} ${LIBDIR}/
