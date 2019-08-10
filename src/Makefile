.PHONY: clean all info

HOME_EXERCISE = 
TARGETS := main mainf
SOURCES := main.cpp mainf.cpp
OBJS    := $(SOURCES:.cpp=.o) $(SOURCES:.cpp=_d.o) $(SOURCES:.cpp=_p.o)


FDDEBUGFLAGS:= -g -O3 
DEBUGFLAGS:= -g -O3 -DNFDDEBUG
RELEASEFLAGS:= -O3 -DNDEBUG -DNFDDEBUG
PARALLELFLAGS:= -D_GLIBCXX_PARALLEL -fopenmp -pthread
CXXFLAGS := -Wall -std=c++17 $(RELEASEFLAGS)
LDFLAGS  := 
LIBS =   -lstdc++ 
OPENMPLIB =-L/usr/lib -fopenmp -pthread

EXAMPLE_DEPS = Makefile


all: $(TARGETS)

clean:
	/bin/rm -f $(OBJS) $(TARGETS)

info:
	@echo Compiler:	 CXX	  = $(CXX)
	@echo Compile command: COMPILE.cc  = $(COMPILE.cc)
	@echo Link command:    LINK.cc     = $(LINK.cc)


main.o: main.cpp $(EXAMPLE_DEPS)
	$(CXX) -c $(CXXFLAGS) $(RELEASEFLAGS)  -o main.o main.cpp
main: main.o
	$(CXX) -o main $(LDFLAGS) main.o $(LIBS)


test:
	g++ -I$(HOME)/install/include test.cpp -o test -L$(HOME)/install/lib -lnlopt -lm

main: main.cpp
	g++ -I$(HOME)/install/include main.cpp $(RELEASEFLAGS) -o main -L$(HOME)/install/lib -lnlopt -lm -Wl,-rpath=$(HOME)/install/lib

mainfd: main.cpp
	g++ -I$(HOME)/install/include main.cpp $(FDDEBUGFLAGS) -o mainfd -L$(HOME)/install/lib -lnlopt -lm -Wl,-rpath=$(HOME)/install/lib

maind: main.cpp
	g++ -I$(HOME)/install/include main.cpp $(DEBUGFLAGS) -o maind -L$(HOME)/install/lib -lnlopt -lm -Wl,-rpath=$(HOME)/install/lib


solvePDE: solvePDE.cpp
	g++ -I$(HOME)/install/include solvePDE.cpp $(RELEASEFLAGS) -o solvePDE -L$(HOME)/install/lib -lnlopt -lm -Wl,-rpath=$(HOME)/install/lib




