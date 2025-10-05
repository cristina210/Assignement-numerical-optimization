# Optimization Methods: Nelder-Mead & Modified Newton

## Overview
This project implements and tests two optimization methods in MATLAB:  

1. **Nelder-Mead Method** – a derivative-free simplex-based algorithm for unconstrained optimization.  
2. **Modified Newton Method** – a Newton-type method enhanced with backtracking line search.  

The implementations are applied to the **Rosenbrock function** and additional benchmark problems to study their behavior, convergence, and performance.

---

## Implemented Features

### 1. Nelder-Mead Method
- Standard Nelder-Mead simplex procedure.  
- Parameter tuning for reflection, expansion, contraction, and shrinkage coefficients.   
- Tested on varying dimensions and various randomly generated points.  

### 2. Modified Newton Method
- Uses exact or finite-difference gradients and Hessians.  
- Backtracking line search with tuning of parameters
- Supports derivative approximation via finite differences:
  - Absolute increments: \( h = 10^{-k}, k = 2,4,6,8,10,12 \)  
  - Relative increments: \( h_i = 10^{-k} |x_i|, k = 2,4,6,8,10,12 \)

---

## Test Problems
### Objective function: 
- Rosenbrock Function
- Wood Function
- Powell Function
### Dimension of space:
- Nelder-Mead: 2, 10, 25, 50
- Modified Newton: 2, 10^3, 10^4, 10^5
### Initial points:
- ten initial random starting point

---

## Comparison Criteria
For each method and problem, the following metrics are reported:

| Metric | Description |
|--------|-------------|
| Success/Failure | Number of runs that satisfy stopping criterion |
| Iterations | Number of iterations to converge |
| Convergence Rate | Experimental rate of convergence |
| Execution Time | time required for convergence |

Tables and plots are generated to visually compare performance across methods, starting points, and dimensions.

---

