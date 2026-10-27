all: ActinMixedNucleates.so

CXX=icpx
PYTHON3=python3

LDLIBS= -qmkl
OPENMP=-qopenmp
PYBINDINCLUDE:=`$(PYTHON3) -m pybind11 --includes`

CXXFLAGS=-Wall -O3 -std=c++11 -fPIC $(OPENMP) $(LDLIBS)  $(PYBINDINCLUDE)

%.so:%.cpp makefile
	$(CXX) -shared $< -o $@ $(CXXFLAGS)

clean:
	rm -f *.so


