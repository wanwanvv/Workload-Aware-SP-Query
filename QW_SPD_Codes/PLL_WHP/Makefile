CXX = g++
CXXFLAGS = -g -w -Wall -Wextra -O3 -std=c++11 -fopenmp

all: bin bin/pll_whp_run

bin/pll_whp_run: main.cpp
	$(CXX) $(CXXFLAGS) -Isrc -o $@ $^

bin:
	mkdir -p bin 

clean:
	rm -rf bin
