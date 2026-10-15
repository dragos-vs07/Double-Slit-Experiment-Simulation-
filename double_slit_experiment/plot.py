import matplotlib.animation as animation
from matplotlib import pyplot as plt
from matplotlib.colors import PowerNorm
import numpy as np
import subprocess

root_path = "C:/Users/Dragos/Desktop/Double Slit Experiment Simulation/double_slit_experiment"
subprocess.run(["g++", root_path + "/main.cpp", "-o", root_path + "/main.exe"], check=True)
subprocess.run([root_path + "/main.exe"], check=True)

files = [root_path + f"/snapshots/output_{n}.txt" for n in range(0, 1001, 10)]
fig, ax = plt.subplots()

data0 = np.loadtxt(files[0])
vmax = 0.05 * max(np.loadtxt(f).max() for f in files)   # try 0.02 to 0.2

im = ax.imshow(data0, origin="lower", extent=[0, 24, 0, 40],
               cmap="magma",
               norm=PowerNorm(gamma=0.4, vmin=0, vmax=vmax))

def update(k):
    im.set_data(np.loadtxt(files[k])); return [im]
ani = animation.FuncAnimation(fig, update, frames=len(files), interval=200)
plt.show()