# COMP670_Final_Project_Spring_2025_Deliverables
My COMP 670 Project Deliverables
## Steps to Reproduce
 1. For [hyperbolic1D_upwind.cpp](https://github.com/csrc-sdsu/mole/blob/master/examples/cpp/hyperbolic1D_upwind.cpp
), since it's already merged in master, you should be able to run `make` in your local mole repository in the directory `examples/cpp/` and then run `./hyperbolic1D_upwind`
    - If it doesn't work, then maybe you could compile the file individually, e.g., what I would do on my machine is, in the directory `examples/cpp/`:
       1. Put the `Makefile` from this repo into the directory `examples/cpp/`
       2. Then compile  `/opt/homebrew/Cellar/llvm/20.1.2/bin/clang++ -std=c++20 -O3 -fopenmp -DARMA_DONT_USE_WRAPPER -Wall -DARMA_SUPERLU_INCLUDE_DIR=\"/opt/homebrew/Cellar/superlu/7.0.0/include\"  -I.  -I/opt/homebrew/Cellar/armadillo/14.4.1/include -I../../src/cpp -I/opt/homebrew/Cellar/superlu/7.0.0/include -o hyperbolic1D_upwind hyperbolic1D_upwind.cpp -L../../src/cpp -lmole -L/opt/homebrew/Cellar/openblas/0.3.29/lib -L/opt/homebrew/Cellar/armadillo/14.4.1/lib -L/opt/homebrew/Cellar/superlu/7.0.0/lib -lopenblas -lsuperlu -L/opt/homebrew/opt/superlu/lib -L/opt/homebrew/Cellar/armadillo/14.4.1 -Wl,-rpath,/opt/homebrew/Cellar/armadillo/14.4.1 -larmadillo -L/opt/homebrew/Cellar/superlu/7.0.0/lib -Wl,-rpath,/opt/homebrew/Cellar/superlu/7.0.0/lib -lsuperlu`
       3. And run `./hyperbolic1D_upwind`
 2. For `elliptic2D_case2.cpp`, I would compile like the following:
     1. Put the `Makefile` from this repo into the directory `examples/cpp/`
     2. Then compile:
         - Run `/opt/homebrew/Cellar/llvm/20.1.2/bin/clang++ -std=c++20 -O3 -fopenmp -DARMA_DONT_USE_WRAPPER -Wall -DARMA_SUPERLU_INCLUDE_DIR=\"/opt/homebrew/Cellar/superlu/7.0.0/include\"  -I.  -I/opt/homebrew/Cellar/armadillo/14.4.1/include -I../../src/cpp -I/opt/homebrew/Cellar/superlu/7.0.0/include -o elliptic2D_case2 elliptic2D_case2.cpp -L../../src/cpp -lmole -L/opt/homebrew/Cellar/openblas/0.3.29/lib -L/opt/homebrew/Cellar/armadillo/14.4.1/lib -L/opt/homebrew/Cellar/superlu/7.0.0/lib -lopenblas -lsuperlu -L/opt/homebrew/opt/superlu/lib -L/opt/homebrew/Cellar/armadillo/14.4.1 -Wl,-rpath,/opt/homebrew/Cellar/armadillo/14.4.1 -larmadillo -L/opt/homebrew/Cellar/superlu/7.0.0/lib -Wl,-rpath,/opt/homebrew/Cellar/superlu/7.0.0/lib -lsuperlu`
         - Or run: `make`
     3. And run `./elliptic2D_case2`
        
           
