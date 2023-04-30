
VERSION=1.6
TCLINC=-I/usr/include/tcl8.6
TCLLIB=-L/usr/lib64 -ltclstub8.6
PLUGINDIR=${HOME}/lib/vmd/plugins/LINUXAMD64/tcl

# Tested with g++ and clang++
#CXX=g++
#CXX=clang++
CPPFLAGS=-fpic -O3 -Wall ${TCLINC} -DVERSION=\"${VERSION}\" -DUSE_TCL_STUBS

all: qwrap.so.${VERSION} pkgIndex.tcl

qwrap.so.${VERSION}: qwrap.o
	$(CXX) -shared qwrap.o -o qwrap.so.${VERSION} ${TCLLIB}
	ln -s -f qwrap.so.${VERSION} qwrap.so

qwrap.tar.gz: qwrap.cpp Makefile pkgIndex.tcl
	tar czf qwrap${VERSION}.tar.gz qwrap.cpp Makefile pkgIndex.tcl

pkgIndex.tcl: qwrap.so
	tclsh mkindex.tcl

clean:
	rm -f *.o *.so

install: all
	mkdir -p ${PLUGINDIR}/qwrap${VERSION} && \
		rsync -av qwrap.so pkgIndex.tcl ${PLUGINDIR}/qwrap${VERSION}/. 
