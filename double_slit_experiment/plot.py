import matplotlib.animation as animation
from matplotlib import pyplot as plt
import numpy as np

files = [f"double_slit_experiment/snapshots/output_{n}.txt" for n in range(0, 101, 10)]
fig, ax = plt.subplots()
im = ax.imshow(np.loadtxt(files[0]), origin="lower", extent=[0, 24, 0, 40])
def update(k):
    im.set_data(np.loadtxt(files[k])); return [im]
ani = animation.FuncAnimation(fig, update, frames=len(files), interval=200)
plt.show()