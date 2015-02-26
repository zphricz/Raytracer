CXXFLAGS = -std=c++11 -Ofast
LDFLAGS = -lSDL2
OS = $(shell uname -s)
SRC = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp, %.o, $(SRC))
HEADERS = $(patsubst %.cpp, %.h, $(SRC))
DEPS = $(patsubst %.cpp, %.d, $(SRC))
ELFNAME = raytracer

ifeq ($(OS), Darwin)
	CXX = clang++
endif
ifeq ($(OS), Linux)
	LDFLAGS += -lpthread
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
