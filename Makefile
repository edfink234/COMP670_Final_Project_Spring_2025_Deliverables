MOLE = ../../src/cpp

ifndef ARMA
	ARMA = /opt/homebrew/Cellar/armadillo/14.4.1
endif

CXX = /opt/homebrew/Cellar/llvm/20.1.2/bin/clang++ -std=c++20
DEBUG = 1

ifeq ($(DEBUG),1)
CXXFLAGS = -g
else
CXXFLAGS = -O3
endif

CXXFLAGS += -fopenmp -DARMA_DONT_USE_WRAPPER -Wall -Wextra -Wmisleading-indentation -Wcomment
CXXFLAGS += -DARMA_SUPERLU_INCLUDE_DIR=\"$(SUPERLU)/include\"

INCPATH = -I. $(if $(ARMA), -I$(ARMA)/include) -I$(MOLE) -I$(SUPERLU)/include

OPENBLAS = /opt/homebrew/Cellar/openblas/0.3.29
ARMADILLO = /opt/homebrew/Cellar/armadillo/14.4.1
SUPERLU = /opt/homebrew/Cellar/superlu/7.0.0
LIBS = -L$(MOLE) -lmole -L$(OPENBLAS)/lib -L$(ARMADILLO)/lib -L$(SUPERLU)/lib -lopenblas -lsuperlu -L/opt/homebrew/opt/superlu/lib

ifdef ARMA
LIBS += -L$(ARMA) -Wl,-rpath,$(ARMA) -larmadillo
else
ifeq (,$(filter clean,$(MAKECMDGOALS)))
$(warning The path to Armadillo's library was not provided. I will look for it in the standard system directories.)
endif
endif

ifdef EIGEN
CXXFLAGS += -DEIGEN
INCPATH += -I$(EIGEN)
else
ifdef SUPERLU
LIBS += -L$(SUPERLU)/lib -Wl,-rpath,$(SUPERLU)/lib -lsuperlu 
else
LIBS += -lsuperlu 
endif
endif


all: transport1D schrodinger1D elliptic1D elliptic2D parabolic1D convection_diffusion RK2 elliptic2D_case2 hyperbolic1D_upwind

transport1D: transport1D.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o transport1D transport1D.cpp $(LIBS)

schrodinger1D: schrodinger1D.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o schrodinger1D schrodinger1D.cpp $(LIBS)

elliptic1D: elliptic1D.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o elliptic1D elliptic1D.cpp $(LIBS)

elliptic2D: elliptic2D.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o elliptic2D elliptic2D.cpp $(LIBS)
	
parabolic1D: parabolic1D.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o parabolic1D parabolic1D.cpp $(LIBS)

convection_diffusion: convection_diffusion3D.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o convection_diffusion convection_diffusion3D.cpp $(LIBS)

RK2: RK2.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o RK2 RK2.cpp $(LIBS)

elliptic2D_case2: elliptic2D_case2.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o elliptic2D_case2 elliptic2D_case2.cpp $(LIBS)

hyperbolic1D_upwind: hyperbolic1D_upwind.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -o hyperbolic1D_upwind hyperbolic1D_upwind.cpp $(LIBS)

clean:
	rm -f transport1D schrodinger1D parabolic1D elliptic1D elliptic2D convection_diffusion RK2 hyperbolic1D_upwind
