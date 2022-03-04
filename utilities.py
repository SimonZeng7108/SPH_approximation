import numpy as np
#consturct weights matrix to eliminate particles outside radius
def Cutoff(r, c):
    if 0 < r < c:
        y = 1
    else:
        y = 0
    return y

def CalculateDensity(m, kernel, x, y, h, c ):
    'Compute density based on particle distances'
    density = np.sum(m * kernel(x, y, h, c))
    # normalised_density = density / np.sum(kernel(distance, h, c))
    return density

def CalculateDensity_taichi(m, kernel, r, h):
    'Compute density based on particle distances'
    density = np.sum(m * kernel(r, h))
    # normalised_density = density / np.sum(kernel(distance, h, c))
    return density
