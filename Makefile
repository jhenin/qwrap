
VERSION=1.5
TCLINC=-I/usr/include/tcl8.6
TCLLIB=-I/usr/lib64
PLUGINDIR=${HOME}/lib/vmd/plugins/LINUXAMD64/tcl

CPP=g++
CPPFLAGS=-fpic -O3 -Wall -ansi -pedantic -fno-for-scope ${TCLINC} -DVERSION=\"${VERSION}\"

all: qwrap.so pkgIndex.tcl

qwrap.so: qwrap.o
	$(CPP) -shared qwrap.o -o qwrap.so ${TCLLIB}

qwrap.tar.gz: qwrap.cpp Makefile pkgIndex.tcl
	tar czf qwrap${VERSION}.tar.gz qwrap.cpp Makefile pkgIndex.tcl

pkgIndex.tcl: qwrap.so
	tclsh mkindex.tcl

clean:
	rm *.o *.so

install: all
	mkdir -p ${PLUGINDIR}/qwrap${VERSION} && \
		rsync -av qwrap.so pkgIndex.tcl ${PLUGINDIR}/qwrap${VERSION}/. 
