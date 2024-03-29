CC = g++
CXXFLAGS = -Wall -Werror -O3 $(shell pkg-config --cflags jellyfish-2.0) -std=c++11 -Wno-unused-local-typedefs
LDFLAGS := $(LDFLAGS) $(shell pkg-config --libs-only-L jellyfish-2.0) -L$(BOOST_ROOT)/lib 
LDLIBS = -lboost_iostreams -lboost_timer -lboost_chrono -lboost_system -lboost_program_options  $(shell pkg-config --libs-only-l jellyfish-2.0)

jellyfish-join: jellyfish-join.cc
clean:
	rm -f jellyfish-join
