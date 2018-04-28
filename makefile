.PHONY:
all: goap

CXXFLAGS := -std=c++17
goap: main.o goap.o astar.o
	$(CXX) $(CXXFLAGS) -o $@ $^

.PHONY:
clean:
	rm goap *.o
