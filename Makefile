CC=g++
# CFLAGS=-std=c++11
CFLAGS=-std=c++11 -O3 -w
# CFLAGS += $(shell $(CC) -fopenmp -E - < /dev/null > /dev/null 2>&1 && echo "-fopenmp")
# LDFLAGS += $(shell $(CC) -fopenmp -E - < /dev/null > /dev/null 2>&1 && echo "-fopenmp")

.PHONY : clean mfe
objects=bin/laf_mfe_vienna bin/laf_mfe_bl*

all: mfe

mfe: src/linearalifold.cpp src/inside.cpp src/outside.cpp src/backtrack.cpp
	mkdir -p bin
	# compile with vienna energy_model
	$(CC) src/linearalifold.cpp src/inside.cpp src/outside.cpp src/backtrack.cpp src/utils/energy_model.cpp $(CFLAGS) -Dlv -Dis_cube_pruning -Dis_candidate_list -DEM_Vienna -o bin/laf_mfe_vienna $(LDFLAGS)
	# compile with bl* energy model
	$(CC) src/linearalifold.cpp src/inside.cpp src/outside.cpp src/backtrack.cpp src/utils/energy_model.cpp $(CFLAGS) -Dlv -Dis_cube_pruning -Dis_candidate_list -DEM_BL_Star -o bin/laf_mfe_bl $(LDFLAGS)

clean:
	-rm $(objects)
