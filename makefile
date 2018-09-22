CC=g++
CFLAGS=-c -std=c++11 -O3 -Wall -Isrc
LDFLAGS=
HEADER = src/util/CommonType.h
SOURCES=main.cpp \
	src/process/BasicProcess.cpp \
    src/util/CSVRead.cpp \
    src/network_model/Network.cpp \
		src/network_model/Connection.cpp \
		src/network_model/Node.cpp \
		src/util/Common.cpp \
		src/math/LinearAlgebra.cpp

OBJECTS=$(SOURCES:.cpp=.o)

EXECUTABLE=molecular

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -O3 -o $@

%.o: %.cpp $(HEADER)
	$(CC) $(CFLAGS) $< -O3 -o $@

clean:
	rm molecular $(OBJECTS)
