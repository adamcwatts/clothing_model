##Numerical Results

#### Boundary Conditions:
T(0mm,t) = T(10mm, t) = 0,  t>0
 
#### Initial Conditions:
T(x,0s) = 35C,  t=0

### Discretion 
L=10mm
100 grid points -> delta_x = 0.1mmm 
delta_t = 0.001 seconds for stability criterion

### Solution Types
1) Analytical solution using Sin and Exponential terms and fourier coefficients (500 terms) to match initial conditions
2) Forward finite difference method 

![SegmentLocal](numerical_scheme_comparisons.gif "Comparison of Numerical Methods")