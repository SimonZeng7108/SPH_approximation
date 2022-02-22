import numpy as np
def RFunction(x, y):
    z = np.sqrt(x**2 + y**2)
    return z

def d_RFunction(x, y):
    r = np.sqrt(x**2 + y**2)
    dx = x / (r+1e-10)
    dy = y / (r+1e-10)
    return dx, dy

def dd_RFunction(x, y):
    r = np.sqrt(x**2 + y**2)
    z = 1 /(r+1e-10)
    return z

def SineFunction(x,y):
    z = np.sin(np.sqrt(x**2 + y**2))
    return z

def d_SineFunction(x,y):
    r = np.sqrt(x**2 + y**2)
    dx = x*np.cos(r)/(r+1e-10)
    dy = y*np.cos(r)/(r+1e-10)
    return dx, dy

def dd_SineFunction(x,y):
    r = np.sqrt(x**2 + y**2)
    ddz = (np.cos(r) - np.sin(r)*r) / (r+1e-10)
    return ddz

def CosFunction(x, y):
    z = -(np.cos(x)**2 + np.cos(y)**2)**2
    return z

def d_CosFunction(x, y):
    dx = 4*(np.cos(x)**3)*np.sin(x) + 4*np.cos(x)*np.sin(x)*np.cos(y)**2
    dy = 4*(np.cos(y)**3)*np.sin(y) + 2*(np.cos(x)**2)*np.sin(2*y)
    return dx, dy

def dd_CosFunction(x, y):
    ddz =2 * (2 * np.cos(2 * x) + np.cos(4 * x) + np.cos(2 * (x - y)) + 2 * np.cos(2 * y) + np.cos(4 * y) + np.cos(2 * (x + y)))
    return ddz

def FrankeFunction(x,y):
    term1 = 0.75*np.exp(-(0.25*(9*x-2)**2) - 0.25*((9*y-2)**2))
    term2 = 0.75*np.exp(-((9*x+1)**2)/49.0 - 0.1*(9*y+1)**2)
    term3 = 0.5*np.exp(-0.25*(9*x-7)**2 - 0.25*((9*y-3)**2))
    term4 = -0.2*np.exp(-(9*x-4)**2 - (9*y-7)**2)
    return term1 + term2 + term3 + term4

def d_FrankeFunction(x, y):
    term1 = 0.75*np.exp(-(0.25*(9*x-2)**2) - 0.25*((9*y-2)**2))
    term2 = 0.75*np.exp(-((9*x+1)**2)/49.0 - 0.1*(9*y+1)**2)
    term3 = 0.5*np.exp(-0.25*(9*x-7)**2 - 0.25*((9*y-3)**2))
    term4 = -0.2*np.exp(-(9*x-4)**2 - (9*y-7)**2)

    dfx = -2*(9*x - 2)*9/4 * term1 - 2*(9*x + 1)*9/49 * term2 
    -2*(9*x - 7)*9/4 * term3 + 2*(9*x - 4)*9 * term4

    dfy = -2*(9*y - 2)*9/4 * term1 - 2*(9*y + 1)*9/10* term2 
    -2*(9*y - 3)*9/4 * term3 + 2*(9*y - 7)*9 * term4
    
    return dfx, dfy

def dd_FrankeFunction(x, y):
    term1 = 0.75*np.exp(-(0.25*(9*x-2)**2) - 0.25*((9*y-2)**2))
    term2 = 0.75*np.exp(-((9*x+1)**2)/49.0 - 0.1*(9*y+1)**2)
    term3 = 0.5*np.exp(-(9*x-7)**2/4.0 - 0.25*((9*y-3)**2))
    term4 = -0.2*np.exp(-(9*x-4)**2 - (9*y-7)**2)
    d_term1_dx = -2*(9*x - 2)*9/4
    d_term1_dy = -2*(9*y - 2)*9/4
    d_term2_dx = -2*(9*x + 1)*9/49
    d_term2_dy = -2*(9*y + 1)*9/10
    d_term3_dx = -2*(9*x - 7)*9/4
    d_term3_dy = -2*(9*y - 3)*9/4
    d_term4_dx = 2*(9*x - 4)*9
    d_term4_dy = 2*(9*y - 7)*9

    dd_term1 = -2*9*9/4 * term1 + d_term1_dx * d_term1_dx * term1 + -2*9*9/4 * term1 + d_term1_dy * d_term1_dy * term1
    dd_term2 = -2*9*9/49 * term2 + d_term2_dx * d_term2_dx * term2 + 0 + d_term2_dy * d_term2_dy * term2
    dd_term3 = -2*9*9/4 * term3 + d_term3_dx * d_term3_dx * term3 + -2*9*9/4 * term3 + d_term3_dy * d_term3_dy * term3
    dd_term4 = 2*9*9 * term4 + d_term4_dx * d_term4_dx * term4 + 2*9*9 * term4 + d_term4_dy * d_term4_dy * term4
    return dd_term1 + dd_term2 + dd_term3 + dd_term4

if __name__ == '__main__':
    from plot_functions import plot_XYZ_3D as plot3D
    from plot_functions import plot_gradient_2D as plot_grad
    function = np.vectorize(FrankeFunction)
    d_function = np.vectorize(d_FrankeFunction)
    dd_function = np.vectorize(dd_FrankeFunction)
    domain = 1                               #grid range
    particles_per_row = 30                   #grid density
    x = np.linspace(-domain, domain, particles_per_row) 
    y = np.linspace(-domain, domain, particles_per_row)
    x, y = np.meshgrid(x.round(decimals=3), np.flipud(y.round(decimals=3)))
    z = function(x, y)                      #calculate function value
    dx, dy = d_function(x, y)               #calculate gradient of function
    dz = np.sqrt(dx**2 + dy**2)             #calculate magnitude of gradient
    ddz = dd_function(x, y)                 #calculate laplacian of function
    fig1 = plot3D(x, y, z, 1, 'Function', keep = True)
    fig2 = plot_grad(x, y, dx, dy, 2, 'Gradient', keep = True)
    fig3 = plot3D(x, y, dz, 3, 'Magnitude of Gradient', keep = True)
    fig4 = plot3D(x, y, ddz, 4, 'Laplacian', show = True)