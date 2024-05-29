CC=     gcc
CXX=        g++
CFLAGS=     -Wall -fopenmp #-fno-inline-functions -fno-inline-functions-called-once
INCLUDES= -Iropebwt2 -Igfatools
CXXFLAGS=   -Wall -D_GLIBCXX_PARALLEL -fopenmp
LIBS=       -lz


#VPATH = ropebwt2

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDES) $< -o $@

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< -o $@

all: CFLAGS+=-g -O2 -DNDEBUG
all: CXXFLAGS+=-g -O3 -DNDEBUG
all: gindex gindexquery

debug: CFLAGS+=-DDEBUG -g -O0
debug: CXXFLAGS+=-DDEBUG -g -O0
debug: gindex gindexquery

gindex: ropebwt2/rle.o ropebwt2/mrope.o ropebwt2/rope.o ropebwt2/rld0.o gfatools/gfa-io.o  gfatools/gfa-base.o  gfatools/gfa-ed.o  gfatools/kalloc.o  main_index.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)


gindexquery: ropebwt2/rle.o ropebwt2/mrope.o ropebwt2/rope.o ropebwt2/rld0.o gfatools/gfa-io.o  gfatools/gfa-base.o  gfatools/gfa-ed.o  gfatools/kalloc.o  main_query.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

clean:
	rm -rf *.o gindex gindexquery
