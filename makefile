CC=g++
CFLAGS=-c -std=c++11 -O3 -Wall -Isrc
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
	$(CC) $(LDFLAGS) $(OBJECTS) -O3 -o $@

%.o: %.cpp $(HEADER)
	$(CC) $(CFLAGS) $< -O3 -o $@

clean:
	rm raam $(OBJECTS)
