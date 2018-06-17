CC=g++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=main.cpp \
	src/process/BasicProcess.cpp \
    src/util/CSVRead.cpp \
    src/model/Node.cpp \
    src/model/Molecule.cpp
 
OBJECTS=$(SOURCES:.cpp=.o)
 
EXECUTABLE=raam
 
all: $(SOURCES) $(EXECUTABLE)
               
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -O -o $@
 
%.o: %.cpp
	$(CC) $(CFLAGS) $< -O -o $@
 
clean:
	rm raam $(OBJECTS)