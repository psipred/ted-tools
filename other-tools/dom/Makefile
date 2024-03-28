CC=gcc
CFLAGS=-I. -O

dom: dom.o util/wt/util.o util/wt/sort.o util/wt/geom.o util/aa/pdbprot.o
	$(CC) -o dom dom.o util/wt/util.o util/wt/sort.o util/wt/geom.o util/aa/pdbprot.o -lm

clean:
	/bin/rm -f dom dom.o util/wt/util.o util/wt/sort.o util/wt/geom.o util/aa/pdbprot.o
