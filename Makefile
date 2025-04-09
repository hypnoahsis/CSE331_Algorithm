# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -O2

# Source files
SRC_DATA = dataset.cpp
SRC_TEST = test.cpp conv_sort.cpp  # Include conv_sort.cpp here

# Executables
BIN_DATA = data_gen
BIN_TEST = test_exec

# Targets
all: $(BIN_DATA) $(BIN_TEST)

$(BIN_DATA): $(SRC_DATA)
	$(CXX) $(CXXFLAGS) -o $(BIN_DATA) $(SRC_DATA)

$(BIN_TEST): $(SRC_TEST)
	$(CXX) $(CXXFLAGS) -o $(BIN_TEST) $(SRC_TEST)

run-data:
	./$(BIN_DATA)

run-test:
	./$(BIN_TEST)

clean:
	rm -f $(BIN_DATA) $(BIN_TEST)


