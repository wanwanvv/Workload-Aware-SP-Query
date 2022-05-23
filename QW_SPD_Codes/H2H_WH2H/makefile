#MakeFile
CXX = g++
CXXFLAGS = -g -O3 -std=c++17
all: bin bin/run

bin/run: main.cpp graph.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^

bin:
	mkdir -p bin

clean:
	rm -rf bin
