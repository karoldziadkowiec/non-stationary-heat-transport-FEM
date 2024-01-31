# Simulation of the non-stationary heat transport process

Technology used in the project: C++ with visual rendering in ParaView

Software for simulating the non-stationary heat transport process using Finite Element Method (FEN) with convection boundary condition - popular method for numerically solving differential equations arising in engineering and mathematical modeling. The software operates in 2D space. The main goal is to calculate the temperature vector {t1} for a specific time step and visualize it on grids. The temperature vector {t1} contains the calculated temperatures at each grid node for a given time step.

Steps taken:
1. Reading mesh data from a text file.
2. Creating a data structure storing integration points and weights according to Gaussian quadrature.
3. Calculation of the matrix [H] for a finite element.
4. Calculation of the matrix [C] for a finite element.
5. Entering the convection boundary condition:
- Calculation of the [HBC] matrix for a finite element.
- Calculation of the {P} vector for a finite element
6. Aggregation of matrices and vectors.
7. Construction of a system of equations for a specific time step.
8. Generating a .vtk file to visualize the temperature on the grid.

MixGrid mesh:
![1](https://github.com/karoldziadkowiec/non-stationary-heat-transport-FEN/blob/master/photos/1.png)
Data:
- Total simulation time: 500 [𝑠]
- Simulation time step: 50 [𝑠]
- Thermal conductivity: 25 [𝑊/𝑚ꞏ 𝐾]
- Convective heat transfer coefficient: 300 [𝑊𝑚²ꞏ 𝐾]
- Ambient temperature: 1200 [°𝐶]
- Initial temperature: 100 [°𝐶]
- Density: 7800 [𝑘𝑔𝑚³]
- Specific heat capacity: 700 [𝐽𝑘𝑔 · 𝐾]
- Number of nodes: 16

Temperature results for integration scheme = 2:
![2](https://github.com/karoldziadkowiec/non-stationary-heat-transport-FEN/blob/master/photos/2.png)
Temperature results for integration scheme = 3:
![3](https://github.com/karoldziadkowiec/non-stationary-heat-transport-FEN/blob/master/photos/3.png)
Temperature results for integration scheme = 4:
![4](https://github.com/karoldziadkowiec/non-stationary-heat-transport-FEN/blob/master/photos/4.png)

50th second of the process:
![50](https://github.com/karoldziadkowiec/non-stationary-heat-transport-FEN/blob/master/photos/50.png)
100 seconds of the process:
![100](https://github.com/karoldziadkowiec/non-stationary-heat-transport-FEN/blob/master/photos/100.png)
150 seconds of the process:
![150](https://github.com/karoldziadkowiec/non-stationary-heat-transport-FEN/blob/master/photos/150.png)
200 seconds of the process:
![200](https://github.com/karoldziadkowiec/non-stationary-heat-transport-FEN/blob/master/photos/200.png)
500 seconds of the process:
![500](https://github.com/karoldziadkowiec/non-stationary-heat-transport-FEN/blob/master/photos/500.png)
