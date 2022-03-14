 
#  SPH Approximation on test functions
This repo uses SPH approximation on different test functions with different kernel functions.<br/>
Test functions: [[see pdf](https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/Derivatives.pdf)] <br/>
- `R_function` <br/>
- `Sin_function` <br/>
- `Cos_function` <br/>

Kernel functions: <br/>
- `Gaussian` <br/>
- `Shepherd` <br/>
- `CubicSpline` <br/>
- `WendlandQuinticC2` <br/>

Approximation formulations: <br/>
- `Monaghan 2005`[[see pdf](https://interactivecomputergraphics.github.io/SPH-Tutorial/pdf/SPH_Tutorial.pdf)] <br/>
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/SPH_form.png" width="300" height="200"><br/>
- `Gambaruto 2015`[[see pdf](https://www.sciencedirect.com/science/article/pii/S0021999115005665)]<br/>
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/gambaruto_form.png" width="300" height="120"><br/>

## Repository overview
[SPH_Approximation.py](https://github.com/SimonZeng7108/SPH_approximation/blob/main/SPH_Approximation.py): Main code to run SPH approximation on test functions  <br/>
[Gambaruto_2015_formulation.py](https://github.com/SimonZeng7108/SPH_approximation/blob/main/Gambaruto_2015_formulation.py): Gambaruto's formulation for approximation  <br/>
[test_functions.py](https://github.com/SimonZeng7108/SPH_approximation/blob/main/test_functions.py): Stores various of test functions <br/>
[kernel_functions.py](https://github.com/SimonZeng7108/SPH_approximation/blob/main/kernel_functions.py): Stores various of kernel functions <br/>
[plot_functions.py](https://github.com/SimonZeng7108/SPH_approximation/blob/main/plot_functions.py): Stores functions for visualisation <br/>
[utilities.py](https://github.com/SimonZeng7108/SPH_approximation/blob/main/utilities.py): Stores utility functions <br/>

## Requirements 
- `numpy`
- `matplotlib`


# Results
*(left: Analytical; right: Numerical)*
## Test function: R function; Kernel: Cubic Spline; Samples: 900:
**Function approximation**<br/>
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/r_function_analytical.png" width="320" height="240">
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/r_function_Numerical.png" width="320" height="240"><br/>

**Gradient approximation**<br/>
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/r_function_gradient.png" width="320" height="240">
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/r_function_gradient_numerical.png" width="320" height="240"><br/>

**Magnitude of Gradient approximation**<br/>
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/r_function_gradient_magnitude.png" width="320" height="240">
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/r_function_gradient_magnitude_numerical.png" width="320" height="240"><br/>

**Laplacian approximation**<br/>
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/r_function_laplacian.png" width="320" height="240">
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/r_function_laplacian_numerical.png" width="320" height="240"><br/>

## Test function: Sin function; Kernel: Cubic Spline; Samples: 900:
**Function approximation**<br/>
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/sin_function_Analytical.png" width="320" height="240">
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/sin_function_numerical.png" width="320" height="240"><br/>

**Gradient approximation**<br/>
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/sin_function_gradient.png" width="320" height="240">
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/sin_function_gradient_numerical.png" width="320" height="240"><br/>

**Magnitude of Gradient approximation**<br/>
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/sin_function_gradient_magnitude.png" width="320" height="240">
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/sin_function_gradient_magnitude_numerical.png" width="320" height="240"><br/>

**Laplacian approximation**<br/>
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/sin_function_laplacian.png" width="320" height="240">
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/sin_function_laplacian_numerical.png" width="320" height="240"><br/>

## Test function: Cos function; Kernel: Cubic Spline; Samples: 900:
**Function approximation**<br/>
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/cos_function_Analytical.png" width="320" height="240">
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/cos_function_numerical.png" width="320" height="240"><br/>

**Gradient approximation**<br/>
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/cos_function_gradient.png" width="320" height="240">
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/cos_function_gradient_numerical.png" width="320" height="240"><br/>

**Magnitude of Gradient approximation**<br/>
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/cos_function_gradient_magnitude.png" width="320" height="240">
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/cos_function_gradient_magnitude_numerica.png" width="320" height="240"><br/>

**Laplacian approximation**<br/>
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/cos_function_laplacian.png" width="320" height="240">
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/cos_function_laplacian_numerical.png" width="320" height="240"><br/>




