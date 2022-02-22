import numpy as np

def Gaussian(x, y, h, c):
    r = np.sqrt(x**2 + y**2)
    alpha = 1/(np.pi * h**2)
    R = -(r/h)**2
    if 0 < r <= c:
        z = alpha * np.exp(R)
    else:
        z = 0.0
    return z

def d_Gaussian(x, y, h, c):
    r = np.sqrt(x**2 + y**2)
    d_alpha = -1/(np.pi * h**4)
    R = -(r/h)**2
    if 0 < r <= c: 
        dz_dx = d_alpha * 2 * x * np.exp(R)
        dz_dy = d_alpha * 2 * y * np.exp(R)
    else:
        dz_dx = 0.0
        dz_dy = 0.0
    return dz_dx, dz_dy

def dd_Gaussian(x, y, h, c):
    r = np.sqrt(x**2 + y**2)
    R = -(r/h)**2
    numerator = 4*(-h**2 + x**2 + y**2)
    denominator = np.pi*h**6
    if 0 < r <= c:
        ddz = numerator / denominator * np.exp(R)
    else:
        ddz = 0.0
    return ddz

# def dd_Gaussian(x, y, h, c):
#     r = np.sqrt(x**2 + y**2)
#     R = -(r/h)**2
#     alpha = 1 / (np.pi * h**2)
#     if 0 < r <= c:
#         ddz = alpha * 4*(h**2 - r**2) * np.exp(R) / (h**4)
#     else:
#         ddz = 0.0
#     return ddz

def Shepherd(x, y, h, c):
    '''Shepherd interpolant kernel'''
    r = np.sqrt(x**2 + y**2)
    if 0 < r < c:
        w = c / (r+1e-10) - 1
    else:
        w = 0.0
    return w

def d_Shepherd(x, y, h, c):
    r = np.sqrt(x**2 + y**2)
    if 0 < r < c:
        dx = c * x / (r**3+1e-10)  - 1
        dy = c * y / (r**3+1e-10)  - 1
    else:
        dx = 0.0
        dy = 0.0
    return dx, dy

def dd_Shepherd(x, y, h, c):
    r = np.sqrt(x**2 + y**2)
    if 0 < r < c:
        ddz = c / (r**3+1e-10) - 1
    else:
        ddz = 0.0 
    return ddz

def CubicSpline(x, y, h, c):
    alpha = 15/(7*np.pi*h**2)
    r = np.sqrt(x**2 + y**2)
    q = r / h
    if 0 < r < c:
        z = alpha  *(2/3 - q**2 + 1/2*q**3)
    elif h <= r < 2*c:
        z = alpha /6 * (2-q)**3
    else:
        z = 0.0
    return z

def d_CubicSpline(x, y, h, c):
    alpha = 15/(7*np.pi*h**2)
    r = np.sqrt(x**2 + y**2)
    q = r / h
    if 0 < r < c:
        dz_dx = alpha * (-2*q+3/2*q**2) * x / (r * h)
        dz_dy = alpha * (-2*q+3/2*q**2) * y / (r * h)
    elif h <= r < 2*c:
        dz_dx = alpha * (-1/2*(2-q)**2) * x / (r * h)
        dz_dy = alpha * (-1/2*(2-q)**2) * y / (r * h)
    else:
        dz_dx = 0.0
        dz_dy = 0.0
    return dz_dx, dz_dy

def dd_CubicSpline(x, y, h, c):
    alpha = 15/(7*np.pi*h**2)
    r = np.sqrt(x**2 + y**2)
    q = r / h
    if 0 < r < c:
        ddz = alpha * (-2 + 3*q) #/ (h**2) #+ alpha *3/2 * (-2*q+3/2*q**2) * 2 / h
    elif h <= r < 2*c:
        ddz = alpha * (2 - q) #/ (h**2) #+ alpha *3/2 * (-2*q+3/2*q**2) * 2 / h
    else:
        ddz = 0.0
    return ddz

def WendlandQuinticC2(x, y, h, c):
    c = 2*h
    r = np.sqrt(x**2 + y**2)
    q = r / h
    alpha = 7/(4*np.pi*h**2)
    if 0 < r <= c:
        z = alpha*(1-q/2)**4*(2*q+1)
    else:
        z = 0.0
    return z

def d_WendlandQuinticC2(x, y, h, c):
    c = 2*h
    r = np.sqrt(x**2 + y**2)
    q = r / h
    alpha = 7/(4*np.pi*h**2)
    if 0 < r <= c:
        dx = alpha*(-2*(1-q/2)**3*(2*q+1)+2*(1-q/2)**4) * x / (r * h)
        dy = alpha*(-2*(1-q/2)**3*(2*q+1)+2*(1-q/2)**4) * y / (r * h)
    else:
        dx = 0.0
        dy = 0.0
    return dx, dy

def dd_WendlandQuinticC2(x, y, h, c):
    c = 2*h 
    r = np.sqrt(x**2 + y**2)
    q = r / h
    alpha = 7/(4*np.pi*h**2)
    if 0 < r <= c:
        ddz = alpha * (3*(1-q/2)**2*(2*q+1) - 8*(1-q/2)**3) / h **2  + alpha*(-2*(1-q/2)**3*(2*q+1)+2*(1-q/2)**4) * 2 / h
    else:
        ddz = 0.0
    return ddz

# def dd_WendlandQuinticC2(x, y, h, c):
#     c = 2*h 
#     r = np.sqrt(x**2 + y**2)
#     q = r / h
#     alpha = 15/(4*np.pi*h**6)
#     if 0 < r <= c:
#         ddz = alpha * (2*h**4 + 4*h*r**3 - 6*r**4) / (r**3)
#     else:
#         ddz = 0.0
#     return ddz


def dd_viscousity(x, y, h, c):
    r = np.sqrt(x**2 + y**2)
    alpha = 40 / (np.pi * h**5)
    if 0 < r <= c:
        ddz = alpha * (h - r)**2
    else:
        ddz = 0
    return ddz

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
    kernel = np.vectorize(Gaussian)
    d_kernel = np.vectorize(d_Gaussian)
    dd_kernel = np.vectorize(dd_Gaussian)
    domain = 2                               #grid range
    particles_per_row = 50                 #grid density
    h = 10 * domain/particles_per_row
    c = 20 * domain/particles_per_row
    x = np.linspace(-domain, domain, particles_per_row)
    y = np.linspace(-domain, domain, particles_per_row)
    x, y = np.meshgrid(x.round(decimals=3), np.flipud(y.round(decimals=3)))
    z = kernel(x, y, h, c)
    print(z.round(3))
    dz_dx, dz_dy = d_kernel(x, y, h, c)
    dz = np.sqrt(dz_dx**2 + dz_dy**2)
    ddz = dd_kernel(x, y, h, c)
    fig = plot3D(x, y, z, 1, 'Function', keep = True)
    fig = plot3D(x, y, dz, 2, 'Gradient', keep = True)
    fig = plot3D(x, y, ddz, 3, 'Laplacian', show = True)
