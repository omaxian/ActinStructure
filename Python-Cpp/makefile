all: ActinMixedNucleates.so

CXX=g++
PYTHON3=python3

LDFLAGS="-L/opt/homebrew/opt/openblas/lib"
CPPFLAGS="-I/opt/homebrew/opt/openblas/include"
LDLIBS= -lopenblas -llapack
OPENMP=-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib -lomp
PYBINDINCLUDE:=`$(PYTHON3) -m pybind11 --includes` -undefined dynamic_lookup

CXXFLAGS=-Wall -O3 -std=c++11 -fPIC $(OPENMP) $(CPPFLAGS) $(LDFLAGS) $(LDLIBS)  $(PYBINDINCLUDE)

%.so:%.cpp makefile
	$(CXX) -shared $< -o $@ $(CXXFLAGS)

clean:
	rm -f *.so #FortranHeap.o


