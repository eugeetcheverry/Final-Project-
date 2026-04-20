import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys

info = np.loadtxt('dw.csv', delimiter = ',', max_rows = 1)
data = np.loadtxt('dw.csv', delimiter = ',', skiprows = 1)

t = np.linspace(0, info[-1], data.shape[1])
norm = np.zeros((data.shape[0]//3, data.shape[1]))
for i in range(norm.shape[0]):
    norm[i][:] = np.linalg.norm(data[i*3:i*3 + 3], axis = 0)


norm_hist = norm.reshape(-1)
t_hist = np.tile(t, data.shape[0]//3)

plt.hist2d(t_hist, norm_hist, bins=[200, 200], cmap = 'gray', vmax = 100)
plt.show()
