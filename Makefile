CC=g++
CFLAGS=-c -Wall -m64 -DDEBUG -g
LDFLAGS=-fPIC
DIR_SRC = ./src
DIR_BIN = ./bin
SOURCES=$(wildcard $(DIR_SRC)/*.cpp) read_oscillation_v01.cpp
TOOLS = $(wildcard ./tools/*.cpp)
TOOL_EXECUTABLES = $(patsubst ./tools/%.cpp,$(DIR_BIN)/%,$(TOOLS))
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=$(DIR_BIN)/read_oscillation_v01

# Edit this for the location of your root executable
ROOTSYS=/usr
CFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
LDFLAGS += $(shell $(ROOTSYS)/bin/root-config --libs) 

CFLAGS += -std=c++17 -I./inc/ -I$(ROOTSYS)/include/ -I/usr/local/include/Minuit2/

all: $(DIR_BIN) $(EXECUTABLE) $(TOOL_EXECUTABLES)

$(DIR_BIN):
	mkdir -p $(DIR_BIN)

$(EXECUTABLE): $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LDFLAGS) -lMinuit2

$(DIR_BIN)/%: ./tools/%.o
	$(CC) -o $@ $< $(LDFLAGS) -lMinuit2

%.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o $(DIR_SRC)/*.o $(DIR_BIN)/* $(OBJECTS)
	rm -f *.pcm *.d *.so
canv:
	rm -f canv*
rroot:
	rm -f *.root
