CC=g++
CFLAGS=-c -std=c++11 -o -O3 -Wall -Isrc
LDFLAGS=
HEADER = src/util/CommonType.h
SOURCES=main.cpp \
		src/ElasticNetworkModel.cpp \
		src/network_model/Connection.cpp \
		src/network_model/Node.cpp \
		src/network_model/HessianMatrix.cpp \
		src/network_model/ReferencePoint.cpp \
    src/util/CSVRead.cpp \
		src/util/Common.cpp

OBJECTS=$(SOURCES:.cpp=.o)

EXECUTABLE=raam

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

%.o: %.cpp $(HEADER)
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm raam $(OBJECTS)
