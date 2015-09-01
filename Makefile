
VERSION=1.0
TCLINC=/usr/include/tcl
PLUGINDIR=${HOME}/lib/vmd/plugins/LINUXAMD64/tcl

CPP=g++
#CPPFLAGS=-fpic -g -I${TCLINC}
CPPFLAGS=-fpic -O3 -I${TCLINC}

all: qwrap.so qwrap.tar.gz

qwrap.so: qwrap.o
	$(CPP) -shared qwrap.o -o qwrap.so 

qwrap.tar.gz: qwrap.cpp Makefile
	tar czf qwrap${VERSION}.tar.gz qwrap.cpp Makefile pkgIndex.tcl

clean:
	rm *.o *.so *.tar.gz

install:
	mkdir -p ${PLUGINDIR}/qwrap${VERSION} && \
		rsync -av qwrap.so pkgIndex.tcl ${PLUGINDIR}/qwrap${VERSION}/. 
