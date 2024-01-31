# Simulation of the non-stationary heat transport process

Technology used in the project: C++ with visual rendering in ParaView

Software for simulating the non-stationary heat transport process using Finite Element Method (FEN) with convection boundary condition - popular method for numerically solving differential equations arising in engineering and mathematical modeling. The software operates in 2D space. The main goal is to calculate the temperature vector {t1} for a specific time step and visualize it on grids. The temperature vector {t1} contains the calculated temperatures at each grid node for a given time step.

Steps taken:
1. Reading mesh data from a text file.
2. Creating a data structure storing integration points and weights according to Gaussian quadrature.
3. Calculation of the matrix [H] for a finite element.
4. Calculation of the matrix [C] for a finite element.
5. Entering the convection boundary condition:
a. Calculation of the [HBC] matrix for a finite element.
b. Calculation of the {P} vector for a finite element
6. Aggregation of matrices and vectors.
7. Construction of a system of equations for a specific time step.
8. Generating a .vtk file to visualize the temperature on the grid.

