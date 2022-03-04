import numpy as np
def cubic_kernel(r_norm, h):
    res = 0.0
    # value of cubic spline smoothing kernel
    k = 40 / 7 / np.pi
    k /= h ** 2
    q = r_norm / h
    if q <= 1.0:
        if q <= 0.5:
            q2 = q * q
            q3 = q2 * q
            res = k * (6.0 * q3 - 6.0 * q2 + 1)
        else:
            res = k * 2 * (1 - q)**3
    # print('cubic kernel ------', res)
    return res


def cubic_kernel_derivative(x, y, h):
    # derivative of cubic spline smoothing kernel
    k = 40 / 7 / np.pi
    k = 6.0 * k / h ** 2
    r = [x, y]
    r_norm = np.sqrt(r[0]**2 + r[1]**2)
    q = r_norm / h
    res = [0.0, 0.0]
    if r_norm > 1e-5 and q <= 1.0:
        grad_q = r / (r_norm * h)
        if q <= 0.5:
            res = k * q * (3.0 * q - 2.0) * grad_q
        else:
            factor = 1.0 - q
            res = k * (-factor * factor) * grad_q
    # print('cubic kernel derivatives------', res)
    return res[0], res[1]

if __name__ == '__main__':
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    from plot_functions import plot_XYZ_3D as plot3D
    from plot_functions import plot_gradient_2D as plot_grad

    # from kernel_functions import Gaussian, d_Gaussian, dd_Gaussian
    # from kernel_functions import Shepherd, d_Shepherd, dd_Shepherd
    # from kernel_functions import CubicSpline, d_CubicSpline, dd_CubicSpline
    # from kernel_functions import WendlandQuinticC2, d_WendlandQuinticC2, dd_WendlandQuinticC2

    #plot kernel function
    kernel = np.vectorize(cubic_kernel)
    d_kernel = np.vectorize(cubic_kernel_derivative)
    domain = 2                               #grid range
    particles_per_row = 50                 #grid density
    h = 20 * domain/particles_per_row
    x = np.linspace(-domain, domain, particles_per_row)
    y = np.linspace(-domain, domain, particles_per_row)
    x, y = np.meshgrid(x.round(decimals=3), np.flipud(y.round(decimals=3)))
    r_norm = np.sqrt(x**2 + y**2)
    z = kernel(r_norm, h)
    dz_dx, dz_dy = d_kernel(x, y, h)
    dz = np.sqrt(dz_dx**2 + dz_dy**2)
    fig = plot3D(x, y, z, 1, 'Function', keep = True)
    fig = plot3D(x, y, dz, 2, 'Gradient', show = True)

