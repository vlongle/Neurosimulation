import numpy as np

filename = "test.txt"
data = np.loadtxt(filename)
y = data[:, 0]
X = data[:, 1:]
