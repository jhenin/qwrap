
VERSION=1.8
TCLINC=-I/usr/include/tcl8.6
TCLLIB=-L/usr/lib64 -ltclstub8.6
PLUGINDIR=${HOME}/lib/vmd/plugins/LINUXAMD64/tcl
VMD?=vmd

# Tested with g++ and clang++
#CXX=g++
#CXX=clang++
CPPFLAGS=-fpic -O3 -Wall ${TCLINC} -DVERSION=\"${VERSION}\" -DUSE_TCL_STUBS

all: qwrap.so.${VERSION} pkgIndex.tcl

qwrap.so.${VERSION}: qwrap.o
	$(CXX) -shared qwrap.o -o qwrap.so.${VERSION} ${TCLLIB}
	ln -s -f qwrap.so.${VERSION} qwrap.so
	cp qwrap.so.${VERSION} Linux_x86_64
	cp qwrap.so Linux_x86_64

qwrap.tar.gz: qwrap.cpp Makefile pkgIndex.tcl
	tar czf qwrap${VERSION}.tar.gz qwrap.cpp Makefile pkgIndex.tcl

pkgIndex.tcl: qwrap.so.${VERSION}
	tclsh mkindex.tcl

clean:
	rm -f *.o *.so

install: all
	mkdir -p ${PLUGINDIR}/qwrap${VERSION} && \
		rsync -av qwrap.so pkgIndex.tcl ${PLUGINDIR}/qwrap${VERSION}/.

test: all
	cd tests && ${VMD} -dispdev text -e test.vmd -eofexit
