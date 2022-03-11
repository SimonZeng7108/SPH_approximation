import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

from utilities import Cutoff, CalculateDensity
from plot_functions import plot_XYZ_3D as plot3D
from plot_functions import plot_gradient_2D as plot_grad
from plot_functions import plot_XYZ_2D as plot2D


#######import functions#######
from test_functions import RFunction, d_RFunction, dd_RFunction
from test_functions import SineFunction, d_SineFunction, dd_SineFunction
from test_functions import CosFunction, d_CosFunction, dd_CosFunction
from test_functions import FrankeFunction, d_FrankeFunction, dd_FrankeFunction
from kernel_functions import Gaussian, d_Gaussian, dd_Gaussian
from kernel_functions import Shepherd, d_Shepherd, dd_Shepherd
from kernel_functions import CubicSpline, d_CubicSpline, dd_CubicSpline
from kernel_functions import WendlandQuinticC2, d_WendlandQuinticC2, dd_WendlandQuinticC2


function = np.vectorize(SineFunction)                      #change the functions here
d_function = np.vectorize(d_SineFunction)
dd_function = np.vectorize(dd_SineFunction)
kernel = np.vectorize(Gaussian)                #change the gradients here
d_kernel = np.vectorize(d_Gaussian)
dd_kernel = np.vectorize(dd_Gaussian)
#######setup constants#######
particles_per_row_list = [30]             #grid density

for particles_per_row in particles_per_row_list:
    domain = 4                                  #grid range
    
    h = 4 * domain/particles_per_row            #smoothing length
    m_v = 0.8 * (domain/particles_per_row)**2   #particle volume   
    a = 4
    c = a * domain/particles_per_row            #cutting off distance
    x = np.linspace(-domain, domain, particles_per_row) 
    y = np.linspace(-domain, domain, particles_per_row)
    x, y = np.meshgrid(x.round(decimals=3), y.round(decimals=3))
    z = function(x, y)                      #calculate function value
    dx, dy = d_function(x, y)               #calculate gradient of function
    dz = np.sqrt(dx**2 + dy**2)             #calculate magnitude of gradient
    ddz = dd_function(x, y)                 #calculate laplacian of function

    #######Plot the analytical function surface#######
    fig1 = plot3D(x, y, z, 1, 'Analytical' , keep = True)
    #######plot gradient of Function#######
    fig2 = plot_grad(x[a : -a, a: -a], y[a : -a, a: -a], dx[a : -a, a: -a], dy[a : -a, a: -a], 2, 'Analytical Gradient', keep = True)
    #######plot magnitude of Gradient#######
    fig3 = plot3D(x[a : -a, a: -a], y[a : -a, a: -a], dz[a : -a, a: -a], 3, 'Analytical Magnitude of Gradient', keep = True)
    # #######plot laplacian#######
    fig4 = plot3D(x[a : -a, a: -a], y[a : -a, a: -a], ddz[a : -a, a: -a], 4, 'Analytical Laplacian', keep = True)
    #######approximate function and its gradient#######
    #initialised density
    density = np.zeros_like(x)
    for i in range(density.shape[0]):
        for j in range(density.shape[1]):
            x_dis = x[i, j] - x
            y_dis = y[i, j] - y
            r = np.sqrt(x_dis**2 + y_dis**2)
            weights = kernel(x_dis, y_dis, h, c)
            density[i, j] = CalculateDensity(m_v, kernel, x_dis, y_dis, h, c)
    mass = np.ones_like(density) * m_v


    cut = np.vectorize(Cutoff)
    approxi_z = np.zeros_like(z)
    approxi_dz_dx = np.zeros_like(dx)
    approxi_dz_dy = np.zeros_like(dy)
    approxi_ddz = np.zeros_like(ddz)
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            x_dis = x[i, j] - x
            y_dis = y[i, j] - y
            r = np.sqrt(x_dis**2 + y_dis**2)

            neighbours = cut(r, c)  
            
            density_i = density[i, j]
            density_j = density * neighbours
            density_j[density_j == 0] = 1e-10
            mass_j = mass * neighbours
            
            x_dis = x_dis * neighbours
            y_dis = y_dis * neighbours
            weights = kernel(x_dis, y_dis, h, c)
            d_weights_dx, d_weights_dy = d_kernel(x_dis, y_dis, h, c)
            dd_weights = dd_kernel(x_dis, y_dis, h, c)
            
            
            phi_i = z[i, j]
            phi_j = z * neighbours
            r[r == 0] = 1e-10

            
            approxi_z[i, j] = np.sum(z * mass_j / density_j * weights)
            approxi_dz_dx[i, j] = density_i*np.sum(mass_j * (phi_i/density_i**2 +  phi_j/density_j**2) * d_weights_dx)
            approxi_dz_dy[i, j] = density_i*np.sum(mass_j * (phi_i/density_i**2 +  phi_j/density_j**2)  * d_weights_dy)
            approxi_ddz[i, j] = 8*np.sum( mass_j / density_j * (phi_j -  phi_i) * r / (r**2) * (np.sqrt(d_weights_dx**2 + d_weights_dy**2)))





    #######plot approximated function values#######
    fig5 = plot3D(x, y, approxi_z, 5, 'Numerical', keep = True)
    #######plot approximated gradient of Function#######
    fig6 = plot_grad(x[a : -a, a: -a], y[a : -a, a: -a], approxi_dz_dx[a : -a, a: -a], approxi_dz_dy[a : -a, a: -a], 6, 'Approximated Gradient', keep = True)
    #######plot magnitude of approximated gradient of Function#######
    approxi_dz = np.sqrt(approxi_dz_dx**2 + approxi_dz_dy**2)
    fig7 = plot3D(x[a : -a, a: -a], y[a : -a, a: -a], approxi_dz[a : -a, a: -a], 7, 'Approximated manitude of Gradient', keep = True)
    #######plot approximated laplacian values#######
    fig8 = plot3D(x[a : -a, a: -a], y[a : -a, a: -a], approxi_ddz[a : -a, a: -a], 8, 'Approximated laplacian', show = True)

    #######caculate RMS-errors#######
    function_rms_error = np.sum(abs((approxi_z[a : -a, a: -a] - z[a : -a, a: -a])))/(x.shape[0]-a)**2
    d_function_rms_error = np.sum(abs((approxi_dz[a : -a, a: -a] - dz[a : -a, a: -a])))/(x.shape[0]-a)**2
    dd_function_rms_error = np.sum(abs((approxi_ddz[a : -a, a: -a] - ddz[a : -a, a: -a])))/(x.shape[0]-a)**2
    print('number of particles:', str(x.shape[0]**2), 
    '   f rms:', str(function_rms_error), 
    '   d_f rms:', str(d_function_rms_error),
    '   dd_f rms:', str(dd_function_rms_error))

