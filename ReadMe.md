# Numerical Solution of Partial Differential Equations

This project is a **homework** for the Numerical Solution of Partial Differential Equations course at **Peking University**.

This project has only been tested on **Linux**.

## Introduction

Currently includes a lightweight **C language library** for sparse matrix computation based on the **Compressed Sparse Row (CSR)** format.
The library supports basic sparse operations, iterative solvers (Jacobi, Gauss-Seidel, Conjugate Gradient), and includes verification examples for the solvers.

***25/11/01 Update:***  
Add a PDE solving part: Enabling solving **2D-Poisson equations** with **Dirichlet Boundary Condition** and **Neumann Boundary Condition**.
Add **python scripts** for visualizing the output data.

***25/12/03 Update:***  
Enhance the PDE solving part: Enabling solving **2D-Parablic PDEs** using **Explicit** and **Implicit (ADI) methods**.
Rearrange the **python scripts** for visualizing time-series data and save animations.

See this project on github: <https://github.com/ZhijunLi-Alessandro/Numerical-Solution-of-PDEs>

## Features

- Basic implementation of Sparse Matrix
  - Efficient **CSR (Compressed Sparse Row)** storage structure  
  - Basic sparse matrix operations: creation, multiplication, destruction  
  - Iterative solvers:
    - Jacobi iteration  
    - Gauss–Seidel iteration  
    - Conjugate Gradient (CG) method  
  - Basic console and file output for CSR sparse matrix
- 2D Poisson Equation Solver
  - Simple orthogonal grid storage structure, provide index translation and region classification
  - Basic grid operations: creation from user-define function, translation, destruction
  - matrix and RHS assembler for 2D Poisson equation with Dirichlet boundary
  - matrix and RHS assembler for 2D Poisson equation with Neumann boundary
- 2D Parabolic Equation Solver
  - Using same grid structure as 2D Poisson Equation Solver
  - Explicit solver
  - ADI Method
- Python scripts for simple visualizing the output results
- Doxygen-compatible documentation

## Project Structure

```
NPDE/
│
├── include/ # Header files
│ ├── bessel.h # Compute Bessel functions
│ ├── csr.h # CSR matrix definition & operations, iterative solvers
│ ├── grid.h # 2D grid definition & operations
│ ├── parabolic.h # 2D Parabolic matrix and RHS assembler
│ ├── poisson2d.h # 2D poisson matrix and RHS assembler
│ ├── utils.h # console output functions & csv file output function
│ └── vec.h # basic vector operations
│
├── src/ # Source implementation
│ ├── math
│ | └── bessel.c
│ ├── pde
│ | ├── grid.c
│ | └── poisson2d.c
│ ├── sparse
│ | ├── csr.c
│ | ├── utils.c
│ | └── vec.c
| └── CMakeLists.txt
│
├── tests/ # tests files (will not be compiled by default)
│ ├── test_csr_3x3.c # Small 3x3 linear system test
│ ├── test_csr_5x5.c # Larger 5x5 linear system test
│ ├── 2D-Poisson.c # 2D Poisson equation matrix generator test
│ └── CMakeLists.txt
|
├── examples/ # Toy problem solverse
│ ├── Dirichlet.c # Solve a toy circular problem with Dirichlet boundary condition
│ ├── Neumann.c # Solve a toy circular problem with Neumann boundary condition
│ ├── Parabolic_Verify.c # Test the reference solution
│ ├── Parabolic_Explicit.c # Solve a toy Parabolic problem using explicit method
│ ├── Parabolic_ADI.c # Solve a toy Parabolic problem using ADI method
│ └── CMakeLists.txt
|
├── secripts/ # Python scripts for visualizing the output solution
│ ├── visualize_utils.py # script for visualizing functions
│ └── visualize_... # script for visualizing certain data
│
├── bin/ # Executable files of the tests and examples (output)
│
├── docs/ # Doxygen-generated documentation (output)
│
├── build/ # CMake-generated files (output)
│
├── lib/ # Library file generated for sparse matrix library (currently DLL) (output)
│
├── results/ # Examples output files (solutions and graphs)
│
├── CMakeLists.txt # Root build configuration
├── Doxyfile # Doxygen configuration
└── README.md
```

## Build Instructions

### Requirements
- GCC (or Clang)
- CMake ≥ 3.10
- Python ≥ 3.12 (numpy and matplotlib)
- (Optional) Doxygen ≥ 1.90 for documentation

### Build the project

```bash
mkdir build && cd build
cmake ..
make
```

After compilation, you will find the executables in:

```bash
bin\
```

### Run Examples

```bash
# Run Dirichlet solver
./bin/Dirichlet.c

# Run Neumann solver
./bin/Neumann.c
```

## Solving Poisson Problem

### Toy problem set up:

Computational domain:  
![image](/results/Poisson/region.png)

Function set up:

$$ -u_{xx}-u_{yy}=f(x,y),\;\;\;(x,y)\in\Omega $$
$$ f(x,y)=\sin(\pi x)\cos(2\pi y) $$

Exact Solution:

$$ u=\frac{1}{5\pi^2}\sin(\pi x)\cos(2\pi y) $$

Exact solution:  
![image](/results/Poisson/Dirichlet_exact.png)

### Dirichlet Boundary

Output of *Dirichlet* gird setup:  
```bash
# Run
[[4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 6 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 6 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 6 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]]
```

Solution visualization:  
![image](/results/Poisson/Dirichlet_solution.png)

Solution Error visualization:  
![image](/results/Poisson/Dirichlet_error.png)

### Neumann Boundary
Output of *Neumann* gird setup (eliminate several corners):  
```bash
[[4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 6 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 6 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 6 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 1 1 1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 1 1 1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [7 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 5 1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ],
 [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]]
```

Solution visualization:  
![image](/results/Poisson/Neumann_solution.png)

Solution Error visualization:  
![image](/results/Poisson/Neumann_error.png)

**Notice:**  
Under current grid point division, the interpolation of boundary points on the boundary in the Dirichlet boundary condition is accurate in the vast majority of edges. However, the approximation of the boundary normal derivative using the first-order difference in Neumann boundary conditions has errors at almost all boundary positions. This might be one of the reasons why the numerical solution error under the Neumann boundary condition is significantly higher than that under the Dirichlet boundary condition. The error on the boundary conditions can also be seen from the part near the boundary in the above figure.

## Solving Parabolic PDE
## Toy problem set up:

Computational domain:  
![image](/results/Poisson/region.png)

Function set up:

$$  \begin{cases}
      \begin{aligned}
        & u_{t} = \Delta u+\rho(t) \delta(x_0, y_0),\;\;\; &(x,y)\in\Omega, \;t>0 \\
        & u(x,y,0) = Re\{ -\frac{1}{4}H_0^{(2)}(\frac{1-i}{\sqrt{2}} \sqrt{(x-x_0)^2+(y-y_0)^2}) \},\;\;\; &(x,y)\in \bar{\Omega}\\
        & u(x,y,t) = Re\{ -\frac{1}{4} e^{it} H_0^{(2)}(\frac{1-i}{\sqrt{2}} \sqrt{(x-x_0)^2+(y-y_0)^2}) \},\;\;\; &(x,y) \in \partial \Omega,\;t>0
      \end{aligned}
    \end{cases}$$
$$  $$

Exact Solution:

$$ u(x,y,t) = Re\{ -\frac{1}{4} e^{it} H_0^{(2)}(\frac{1-i}{\sqrt{2}} \sqrt{(x-x_0)^2+(y-y_0)^2}) \} $$

Numerical solution (Approximate Bessel functions using series):  
![image](/results/Parabolic/solution_refined.gif)

### Explicit method

```bash
# run case
bin/Parabolic_Explicit
# visualize the solution
python scripts/visualize_Parabolic_Explicit.py
```

Solution visualization:  
![image](/results/Parabolic/Explicit.gif)

Solution Error (Removed the cneter point $(x_0,y_0)$) visualization:  
![image](/results/Parabolic/Explicit_error_without_center.gif)

See also the solution for Refined mesh and its error:  
![image](/results/Parabolic/Explicit_refined.gif)![image](/results/Parabolic/Explicit_error_refined_without_center.gif)

### Implicit method

```bash
# run case
bin/Parabolic_ADI
# visualize the solution
python scripts/visualize_Parabolic_ADI.py
```

Solution visualization:  
![image](/results/Parabolic/ADI.gif)

Solution Error (Removed the cneter point $(x_0,y_0)$) visualization:  
![image](/results/Parabolic/ADI_error_without_center.gif)

See also the solution for Refined mesh and its error:  
![image](/results/Parabolic/ADI_refined.gif)![image](/results/Parabolic/ADI_error_refined_without_center.gif)

## Doxygen Documentation
 
Current version doxygen documents are contained in `docs\`.

To generate documentation:

```bash
doxygen Doxyfile
```

Then Open:

```
docs/html/index.html
```

## Future Extenions

- Compatible with boundary condition
- Other iterative solvers

## Autor

**李智俊** ***(Zhijun Li)***  
Email: [<zhijunli25@stu.pku.edu.cn>]  
Project start: 2025-09  
License: MIT  