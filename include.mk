#Modify this variable to set the location of sonLib
sonLibRootDir?=${rootPath}/../sonLib
sonLibDir=${sonLibRootDir}/lib
#Use sonLib bin and lib dirs
BINDIR ?= ${sonLibRootDir}/bin
LIBDIR ?= ${sonLibDir}
CPPFLAGS +=-I${sonLibRootDir}/C/inc -I${sonLibRootDir}/externalTools/cutest/

# Hack to include openmp on os x after "brew install lomp
CFLAGS+= -Xpreprocessor -fopenmp -lomp

include  ${sonLibRootDir}/include.mk

LDLIBS = ${sonLibDir}/sonLib.a ${sonLibDir}/cuTest.a ${dblibs} ${LIBS}
LIBDEPENDS = ${sonLibDir}/sonLib.a ${sonLibDir}/cuTest.a 

