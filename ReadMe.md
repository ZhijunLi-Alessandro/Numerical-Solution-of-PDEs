# Numerical Solution of Partial Differential Equations

This project is a **homework** for the Numerical Solution of Partial Differential Equations course at **Peking University**.

This project has only been tested on **Linux**.

## Introduction

Currently includes a lightweight **C language library** for sparse matrix computation based on the **Compressed Sparse Row (CSR)** format.
The library supports basic sparse operations, iterative solvers (Jacobi, Gauss-Seidel, Conjugate Gradient), and includes verification examples for the solvers.

The current project includes an example of sparse matrix generation for solving the Poisson equation within a regular rectangular boundary, which will be updated later for the complete problem-solving.

See this project on github: <https://github.com/ZhijunLi-Alessandro/Numerical-Solution-of-PDEs>

## Features

- Efficient **CSR (Compressed Sparse Row)** storage structure  
- Basic sparse matrix operations: creation, multiplication, destruction  
- Iterative solvers:
  - Jacobi iteration  
  - Gauss–Seidel iteration  
  - Conjugate Gradient (CG) method  
- Basic console output for CSR sparse matrix
- Doxygen-compatible documentation and examples  

## Project Structure

```
NPDE/
│
├── include/ # Header files
│ ├── csr.h # CSR matrix definition & operations, iterative solvers
│ ├── utils.h # console output functions
│ └── vec.h # basic vector operations
│
├── src/ # Source implementation
│ ├── csr.c
│ ├── utils.c
│ ├── vec.c
| └── CMakeLists.txt
│
├── tests/ # Validation and example files
│ ├── test_csr_3x3.c # Small 3x3 linear system test
│ ├── test_csr_5x5.c # Larger 5x5 linear system test
│ └── CMakeLists.txt
|
├── examples/ # Temprary 2D-Poisson equaion sparse matrix generater
│ ├── 2D-Poisson.c # 2D Poisson equation matrix generator test
│ └── CMakeLists.txt
│
├── bin/ # Executable files of the tests and examples (output)
│
├── docs/ # Doxygen-generated documentation (output)
│
├── build/ # CMake-generated files (output)
│
├── lib/ # Library file generated for sparse matrix library (currently DLL) (output)
│
├── CMakeLists.txt # Root build configuration
├── Doxyfile # Doxygen configuration
└── README.md
```

## Build Instructions

### Requirements
- GCC (or Clang)
- CMake ≥ 3.10
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
# Run 3x3 CSR test
./bin/test_csr_3x3

# Run 5x5 CSR test
./bin/test_csr_5x5

# Run 2D Poisson test, <nx> points in x direction and <ny> poinst in y direction
./bin/2D-Poisson <nx> <ny>
```

### Example Output

Output of *iterative solvers*:

```bash
# Run
./bin/test_csr_5x5

# Output
Matrix A in CSR format:
[[2.000000 -1.000000 0.000000 0.000000 0.000000 ],
 [-1.000000 2.000000 -1.000000 0.000000 0.000000 ],
 [0.000000 -1.000000 2.000000 -1.000000 0.000000 ],
 [0.000000 0.000000 -1.000000 2.000000 -1.000000 ],
 [0.000000 0.000000 0.000000 -1.000000 2.000000 ]]
Right-hand side vector b:
[1.000000 2.000000 3.000000 4.000000 5.000000 ]
Jacobi Solution x:
[5.828066 10.657636 13.489464 13.324303 9.161399 ]
Gauss-Seidel Solution x:
[5.833327 10.666657 13.499991 13.333326 9.166663 ]
Conjugate Gradient Solution x:
[5.833333 10.666667 13.500000 13.333333 9.166667 ]
```

Output of *matrix generator for 2D Poisson equation*:

```bash
# Run
./bin/2D-Poisson 3 3

# Output
Generated 2D Poisson matrix in CSR format:
[[1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 -1.0 0.0 0.0 -1.0 4.0 -1.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 -1.0 0.0 0.0 -1.0 4.0 -1.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 -1.0 4.0 -1.0 0.0 0.0 -1.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 -1.0 4.0 -1.0 0.0 0.0 -1.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 ]]
Diagonal matrix D:
[[1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 ]]
Lower triangular matrix L:
[[0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 -1.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 -1.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ]]
Upper triangular matrix U:
[[0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 -1.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 -1.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -1.0 0.0 0.0 -1.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ],
 [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 ]]
```

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

- Numerical solver for Poisson equations
  - Handle non-rectangular boundaries
  - Handle three types of boundary conditions
- possibly add other iterative solvers

## Autor

**李智俊** ***(Zhijun Li)***  
Email: [<zhijunli25@stu.pku.edu.cn>]  
Project start: 2025-09  
License: MIT  