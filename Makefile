cc := g++
all=splitBcBin
DFLAG=-I/hdf5_1_10_5/include /hdf5_1_10_5/lib/libhdf5.so
$(all):
	$(cc) -fopenmp -O3 -std=c++11 main.cpp -o $@ $(DFLAG)
clean:
	rm  splitBcBin splitBcBin_debug
gdb:
	$(cc) -fopenmp -g -std=c++11 main.cpp -o splitBcBin_debug $(DFLAG)