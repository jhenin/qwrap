
VERSION=1.3
TCLINC=/usr/include/tcl8.6
PLUGINDIR=${HOME}/lib/vmd/plugins/LINUXAMD64/tcl

CPP=g++
#CPPFLAGS=-fpic -g -I${TCLINC} -DVERSION=\"${VERSION}\"
CPPFLAGS=-fpic -O3 -I${TCLINC} -DVERSION=\"${VERSION}\"

all: qwrap.so pkgIndex.tcl

qwrap.so: qwrap.o
	$(CPP) -shared qwrap.o -o qwrap.so 

qwrap.tar.gz: qwrap.cpp Makefile pkgIndex.tcl
	tar czf qwrap${VERSION}.tar.gz qwrap.cpp Makefile pkgIndex.tcl

pkgIndex.tcl: qwrap.so
	tclsh mkindex.tcl

clean:
	rm *.o *.so *.tar.gz

install: all
	mkdir -p ${PLUGINDIR}/qwrap${VERSION} && \
		rsync -av qwrap.so pkgIndex.tcl ${PLUGINDIR}/qwrap${VERSION}/. 
