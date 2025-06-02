CXX = g++ -m64
CXXFLAGS = -std=c++17 -Wall -Wextra -O3 -g

OBJS = main.o tsp_parser.o christofides.o held_karp.o \
       novel_algorithm_utils.o novel_algorithm.o

EXEC = tsp_solver

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(OBJS) -lm

main.o: main.cpp tsp_parser.h christofides.h held_karp.h novel_algorithm.h novel_algorithm_utils.h
	$(CXX) $(CXXFLAGS) -c main.cpp -o main.o

tsp_parser.o: tsp_parser.cpp tsp_parser.h
	$(CXX) $(CXXFLAGS) -c tsp_parser.cpp -o tsp_parser.o

christofides.o: christofides.cpp christofides.h tsp_parser.h
	$(CXX) $(CXXFLAGS) -c christofides.cpp -o christofides.o

held_karp.o: held_karp.cpp held_karp.h tsp_parser.h
	$(CXX) $(CXXFLAGS) -c held_karp.cpp -o held_karp.o

novel_algorithm_utils.o: novel_algorithm_utils.cpp novel_algorithm_utils.h tsp_parser.h
	$(CXX) $(CXXFLAGS) -c novel_algorithm_utils.cpp -o novel_algorithm_utils.o

novel_algorithm.o: novel_algorithm.cpp novel_algorithm.h novel_algorithm_utils.h tsp_parser.h
	$(CXX) $(CXXFLAGS) -c novel_algorithm.cpp -o novel_algorithm.o

clean:
	rm -f $(OBJS) $(EXEC)

.PHONY: all clean

