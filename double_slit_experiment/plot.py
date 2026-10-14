import matplotlib.animation as animation
from matplotlib import pyplot as plt
import numpy as np

file = "double_slit_experiment\output.txt"
fig, ax = plt.subplots()
im = ax.imshow(np.loadtxt(file), origin="lower", extent=[0, 24, 0, 40])
def update(k):
    im.set_data(np.loadtxt(file)); return [im]
ani = animation.FuncAnimation(fig, update, frames=1, interval=200)
plt.show()