## <center> Code and Data for paper "An aggregation-based two-grid method for multilevel block Toeplitz linear system''

Here is the data and MATLAB code for the numerical experiments in paper named "An aggregation-based two-grid method for multilevel block Toeplitz linear system''. The numerical experiments contain two parts: the artifical examples and the practical cases. 

---

### 1. Preparations
Before run the experiments, you should check two things:
1. Please confirm that your MATLAB includes the [Global Optimization Toolbox](https://www.mathworks.com/products/global-optimization.html).
2. Add the "help function" and "data" folder to the MATLAB working directory.

---

### 2. Artifical examples
Run the file [test_artifical.m](test_artifical.m) to get the numerical results for the artifical examples which contain the well-conditioned cases and ill-conditioned cases.

---
### 3. Artifical examples
Run the file [test_real.m](test_real.m) to get the numerical results for the practical cases which contain: 
1. The $Q_s$ Lagrangian FEM stiffness matrices for the second order elliptic differential problem on $[0,1]$;
2. The system from the discretization by staggered DG methods of the incompressible Navier-Stokes equations;
3. The system produced in the Flat Panel Display simulation process.

After running the MATLAB code, you can get the iteration steps and running time for each cases. Also, you can change the variable size_test and case_num to test different case with different size. 