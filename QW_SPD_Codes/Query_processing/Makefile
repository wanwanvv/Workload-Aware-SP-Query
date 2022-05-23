CXX = g++
CXXFLAGS = -g -w -Wall -Wextra -std=c++11 -fopenmp

all: bin bin/query_run

bin/query_run: main.cpp
	$(CXX) $(CXXFLAGS) -Isrc -o $@ $^

bin:
	mkdir -p bin 

clean:
	rm -rf bin
