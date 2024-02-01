CXXFLAGS = -O3

all: LinAliFold CentroidLinAliFold

SOURCES1=$(wildcard src/LinAliFold/*.cpp)

LinAliFold: $(SOURCES1)

	$(CXX) $(CXXFLAGS) -o ./bin/LinAliFold $(SOURCES1) -std=c++11 

SOURCES2=$(wildcard src/CentroidLinAliFold/*.cpp)

CentroidLinAliFold: $(SOURCES2) 

	$(CXX) $(CXXFLAGS) -o ./bin/CentroidLinAliFold $(SOURCES2) -std=c++11 