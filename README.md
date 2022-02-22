 
#  SPH Approximation on test functions
This repo uses SPH approximation on different test functions with different kernel functions.<br/>
This is an extended work and Pytorch implementation of [Coral Density Estimation Project](https://vilab.blogs.bristol.ac.uk/2021/01/analysis-of-coral-using-deep-learning/).<br/>
The dataset is downloaded from the original repo by [Ainsley Rutterford](https://github.com/ainsleyrutterford/deep-learning-coral-analysis/blob/master/README.md). 

## Repository overview
[SPH_Approximation.py](https://github.com/SimonZeng7108/SPH_approximation/blob/main/SPH_Approximation.py): Main code to run SPH approximation on test functions <br/>
[Gambaruto_2015_formulation.py](https://github.com/SimonZeng7108/SPH_approximation/blob/main/Gambaruto_2015_formulation.py): Gambaruto's formulation for approximation <br/>
[test_functions.py](https://github.com/SimonZeng7108/SPH_approximation/blob/main/test_functions.py): Stores various of test functions <br/>
[kernel_functions.py](https://github.com/SimonZeng7108/SPH_approximation/blob/main/kernel_functions.py): Stores various of kernel functions <br/>
[plot_functions.py](https://github.com/SimonZeng7108/SPH_approximation/blob/main/plot_functions.py): Stores functions for visualisation <br/>
[utilities.py](https://github.com/SimonZeng7108/SPH_approximation/blob/main/utilities.py): Stores functions for visualisation <br/>

## Requirements 
- `numpy`
- `matplotlib`


# Results
## Test function: R function; Kernel: Gaussian; Samples: 900:
Gradient approximation<br/>
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/r_function_gradient.png" width="300" height="300">
<img src="https://github.com/SimonZeng7108/SPH_approximation/blob/main/results/r_function_gradient_numerical.png" width="300" height="300"><br/>



