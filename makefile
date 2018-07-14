CC=g++
#CFLAGS=-c -std=c++11 -Wall -Isrc
# DEBUG
CFLAGS=-g -c -std=c++11 -Wall -Isrc
LDFLAGS=
HEADER = src/util/Common.h
SOURCES=main.cpp \
	src/process/BasicProcess.cpp \
    src/util/CSVRead.cpp \
    src/network_model/Network.cpp \
		src/network_model/Connection.cpp \
		src/network_model/Node.cpp

OBJECTS=$(SOURCES:.cpp=.o)

EXECUTABLE=molecular

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -O -o $@

%.o: %.cpp $(HEADER)
	$(CC) $(CFLAGS) $< -O -o $@

clean:
	rm molecular $(OBJECTS)
