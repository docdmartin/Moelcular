CC=clang
OPTIMIZATION=-O3
CFLAGS=-c -std=c++11 -Wall -Wextra -Werror -Wredundant-decls -Isrc
LDFLAGS= -L/usr/lib -lstdc++ -lm -pthread
SOURCES=main.cpp \
    Atom.cpp \
		CSVRead.cpp \
		Solver.cpp

OBJECTS=$(SOURCES:.cpp=.o)

EXECUTABLE=potential_response

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(OPTIMIZATION) -o $@

%.o: %.cpp $(HEADER)
	$(CC) $(CFLAGS) $< $(OPTIMIZATION) -o $@

debug: CFLAGS += -DDEBUG -g -Wcast-align -Wmissing-declarations -Wredundant-decls -Wswitch-default -Winvalid-pch -Wformat=2 -Wmissing-format-attribute -Wformat-nonliteral #-Werror=misleading-indentation
debug: OPTIMIZATION = -g
debug: CC = g++
debug: all

release: CFLAGS += -O3
release: OPTIMIZATION = -O3
release: CC = g++
release: all

clean:
	rm $(EXECUTABLE) $(OBJECTS)
