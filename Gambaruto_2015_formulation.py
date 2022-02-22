import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from plot_functions import plot_XYZ_3D as plot3D
from plot_functions import plot_gradient_2D as plot_grad

from test_functions import RFunction, d_RFunction, dd_RFunction
from test_functions import SineFunction, d_SineFunction, dd_SineFunction
from test_functions import CosFunction, d_CosFunction, dd_CosFunction
from test_functions import FrankeFunction, d_FrankeFunction
from kernel_functions import Gaussian, d_Gaussian, dd_Gaussian
from kernel_functions import Shepherd, d_Shepherd, dd_Shepherd
from kernel_functions import CubicSpline, d_CubicSpline, dd_CubicSpline
from kernel_functions import WendlandQuinticC2, d_WendlandQuinticC2, dd_WendlandQuinticC2


#####setup constants#####
testfunction = np.vectorize(CosFunction)      #test function
d_testfunction = np.vectorize(d_CosFunction)  #gradient of test function
dd_function = np.vectorize(dd_CosFunction)    #laplacian of the test function

particles_per_row_list = [40]#[10, 14, 20, 24, 30, 34, 40, 44, 50]                  #grid density
total_number = [i**2 for i in particles_per_row_list]
kernel_list = [Gaussian] #[Gaussian, Shepherd, CubicSpline, WendlandQuinticC2]
Gradient_rms_log = np.zeros((1, len(particles_per_row_list)))
Laplacian_rms_log = np.zeros((1, len(particles_per_row_list)))


for k in kernel_list:
    Gradient_rms = np.array([])
    Laplacian_rms = np.array([])
    for particles_per_row in particles_per_row_list:
        kernel = np.vectorize(k)
        domain = 4                                  #grid range
        h = 4 * domain / particles_per_row
        c = 6 * domain / particles_per_row            #cutting off distance
        x = np.linspace(-domain, domain, particles_per_row) 
        y = np.linspace(-domain, domain, particles_per_row)
        x, y = np.meshgrid(x.round(decimals=3), y.round(decimals=3))
        z = testfunction(x, y)                      #calculate test function value in domain
        dx, dy = d_testfunction(x, y)               #calculate gradient of function
        dz = np.sqrt(dx**2 + dy**2)                 #calculate magnitude of gradient
        ddz = dd_function(x, y)                     #calculate laplacian of function

        #######Plot the analytical function surface#######
        fig1 = plot3D(x, y, z, 1, 'Analytical test function', keep = False)
        #######plot gradient of Function#######
        fig2 = plot_grad(x, y, dx, dy, 2, 'Analytical Gradient', keep = False)
        #######plot magnitude of Gradient#######
        fig3 = plot3D(x, y, dz, 3, 'Analytical Magnitude of Gradient', keep = False)
        #######plot magnitude of Laplacian#######
        fig4 = plot3D(x, y, ddz, 4, 'Analytical Laplacian', keep = True)

        #######approximate function and its gradient#######
        def cutoff(r, c):
            'compute the neighbours weights matrix'
            if r < c:
                y = 1
            else:
                y = 0
            return y
        cut = np.vectorize(cutoff)

        approxi_dz_dx = np.zeros_like(dx)
        approxi_dz_dy = np.zeros_like(dy)
        approxi_ddz = np.zeros_like(ddz)
        for i in range(x.shape[0]):
            for j in range(x.shape[1]):
                x_dis = x[i, j] - x                 #x distance matrix to point of interest
                y_dis = y[i, j] - y                 #y distance matrix to point of interest
                r = np.sqrt(x_dis**2 + y_dis**2)    #norm 2 distance matrix beween point of interest and all other points in the domain
                neighbours = cut(r, c)              #neightbourhood matrix based on r 
                kernel_weights = kernel(x_dis, y_dis, h, c)
                phi_i = z[i, j]                     #scalar value of point of interests
                phi_j = z * neighbours              #scalar value matrix of neighbour points
                x_vector = x_dis * neighbours       #x vector matrix of neighour points
                y_vector = y_dis * neighbours       #y vector matrix of neighour points
                r_squared = r**2 * neighbours
                r_squared[r_squared == 0] = 1e-10

                approxi_dz_dx[i, j] = -2 / np.sum(kernel_weights) * np.sum(kernel_weights*(phi_j-phi_i)*x_vector/r_squared)  #x component of gradient
                approxi_dz_dy[i, j] = -2 / np.sum(kernel_weights) * np.sum(kernel_weights*(phi_j-phi_i)*y_vector/r_squared)  #y component of gradient

                approxi_ddz[i, j] =  4 / np.sum(kernel_weights*r_squared) * np.sum(kernel_weights*(phi_j-phi_i))


        approxi_dz = np.sqrt(approxi_dz_dx**2 + approxi_dz_dy**2)


        #######plot gradient of Function#######
        fig5 = plot_grad(x, y , approxi_dz_dx, approxi_dz_dy, 5, 'Approximated Gradient', keep = True)
        #######plot magnitude of Gradient#######
        fig6 = plot3D(x, y, approxi_dz, 6, 'Approximated Magnitude of Gradient', keep = True)
        #######plot magnitude of Laplacian#######
        fig7 = plot3D(x, y, approxi_ddz, 7, 'Approximated Laplacian', show = True)


        #######caculate RMS-errors#######
        d_function_rms_error = np.sum(abs((approxi_dz - dz)))/(x.shape[0] - 4)
        dd_function_rms_error = np.sum(abs((approxi_ddz - ddz)))/(x.shape[0] - 4)
        print(
        'kernel:', str(k),
        '   number of particles:', str(x.shape[0]**2),  
        '   gradient rms:', str(d_function_rms_error),
        '   laplacian rms:', str(dd_function_rms_error),
        )

        Gradient_rms = np.append(Gradient_rms, d_function_rms_error)
        Laplacian_rms = np.append(Laplacian_rms, dd_function_rms_error)
    Gradient_rms = Gradient_rms.reshape((1, len(Gradient_rms)))
    Laplacian_rms = Laplacian_rms.reshape((1, len(Laplacian_rms)))

    Gradient_rms_log = np.append(Gradient_rms_log, Gradient_rms, axis = 0)
    Laplacian_rms_log = np.append(Laplacian_rms_log, Laplacian_rms, axis = 0)
    
    Gradient_rms = np.array([])
    Laplacian_rms = np.array([])


plt.figure(1)
plt.plot(total_number, Gradient_rms_log[1, :], label = 'Gaussian', color = 'green')
plt.plot(total_number, Gradient_rms_log[2, :], label = 'Shepherd', color = 'orange')
plt.plot(total_number, Gradient_rms_log[3, :], label = 'CubicSpline', color = 'red')
plt.plot(total_number, Gradient_rms_log[4, :], label = 'WendlandQuinticC2', color = 'blue')
plt.legend(loc="upper right")
plt.title('Gradient of Cos Function Approximation')
plt.xlabel("Samples")
plt.ylabel("RMS-error")

plt.figure(2)
plt.plot(total_number, Laplacian_rms_log[1, :], label = 'Gaussian', color = 'green')
plt.plot(total_number, Laplacian_rms_log[2, :], label = 'Shepherd', color = 'orange')
plt.plot(total_number, Laplacian_rms_log[3, :], label = 'CubicSpline', color = 'red')
plt.plot(total_number, Laplacian_rms_log[4, :], label = 'WendlandQuinticC2', color = 'blue')
plt.legend(loc="upper right")
plt.title('Laplacian of Cos Function Approximation')
plt.xlabel("Samples")
plt.ylabel("RMS-error")
plt.show()