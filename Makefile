CXXFLAGS = -std=c++11 -Ofast -Wall -Werror
LDFLAGS = -lSDL2
OS = $(shell uname -s)
SRC = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp, %.o, $(SRC))
HEADERS = $(patsubst %.cpp, %.h, $(SRC))
DEPS = $(patsubst %.cpp, %.d, $(SRC))
ELFNAME = raytracer

ifeq ($(OS), Darwin)
	CXX = g++-4.9
endif
ifeq ($(OS), Linux)
	LDFLAGS += -lpthread
	LDFLAGS += -latomic
	CXX = g++
endif

all: $(ELFNAME)

$(ELFNAME): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o$@ $^ $(LDFLAGS) 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -MMD -MP $< -o $@

-include $(DEPS)

clean:
	rm -f *.d *.o $(ELFNAME)
