# Build and run your ColorFull programs.
#
#   make run EX=01_basic_color_factors   compile examples/<EX>.cc and run it
#   make EX=01_basic_color_factors       just compile it (binary goes in bin/)
#   make clean                           remove compiled programs

EX      ?= 01_basic_color_factors
PREFIX  := $(CURDIR)/install
CXX     := clang++
CXXFLAGS:= -std=c++11 -O2 -I$(PREFIX)/include
LDFLAGS := -L$(PREFIX)/lib/ColorFull -lColorFull -Wl,-rpath,$(PREFIX)/lib/ColorFull

bin/$(EX): examples/$(EX).cc
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

run: bin/$(EX)
	@mkdir -p results && cd results && ../bin/$(EX)

clean:
	rm -rf bin

.PHONY: run clean
