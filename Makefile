CC=g++
CC+=-DDEBUG -g  
CFLAGS=-c -Wall -m64
LDFLAGS=-fPIC
DIR_SRC = ./src
SOURCES=read_oscillation_v01.cxx $(wildcard $(DIR_SRC)/*.cxx)
TOOLS = $(wildcard ./tools/*.cxx)
OBJECTS=$(SOURCES:.cxx=.o)
EXECUTABLE=read_oscillation_v01

# Edit this for the location of your root executable
ROOTSYS=/opt/root/6.22.02-install
CFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
LDFLAGS += $(shell $(ROOTSYS)/bin/root-config --libs) 

CFLAGS += -std=c++17 -I./inc/ -I$(ROOTSYS)/include/ -I/usr/local/include/Minuit2/

all: $(SOURCES) $(EXECUTABLE) $(TOOLS)

$(EXECUTABLE):$(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LDFLAGS) -lMinuit2

#.C.o:
%.o:%.cxx
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o $(DIR_SRC)/*.o; rm $(EXECUTABLE); rm -f *.pcm *.d *.so
canv:
	rm -f canv*
rroot:
	rm -f *.root
