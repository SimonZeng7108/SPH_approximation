from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

def plot_XYZ_3D(X, Y, Z, i, title, show = False, save = False, keep = False):
    '''
    x = x meshgrid, y = y meshgrid, z = z value,
    i = figure number,
    title = figure window name
    '''
    fig = plt.figure(i)
    fig.canvas.set_window_title(title)
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(
        X, 
        Y, 
        Z, 
        cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
    # Customize the z axis.
    # ax.set_zlim(-0.10, 1.40)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    # Add a color bar which maps values to colors.
    fig.colorbar(surf, shrink=0.5, aspect=5)#
    if show == True:
        plt.show()
    if save == True:
        fig.savefig('{}'.format(title))
    if keep == False:
        plt.close()


def plot_gradient_2D(x, y, dx, dy, i, title, show=False, save = False, keep = False):
    fig = plt.figure(i)
    fig.canvas.set_window_title(title)
    plt.quiver(x, y, dx, dy)
    if show == True:
        plt.show()
    if save == True:
        fig.savefig('{}'.format(title))
    if keep == False:
        plt.close()

def plot_XYZ_2D(x, y, z, i, title, show=False, save = False, keep = False):
    fig = plt.figure(i)
    fig.canvas.set_window_title(title)
    plt.scatter(x, y, c=z, cmap='jet',vmin=0, vmax=1, marker = '.')
    plt.colorbar().ax.set_ylabel('z', rotation=270)
    if show == True:
        plt.show()
    if save == True:
        fig.savefig('{}'.format(title))
    if keep == False:
        plt.close()