import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

## Nombre de points
N = 20

## Temps
t = np.linspace(0, 20, 300)

## Création des trajectoires différentes
trajectories = []

for i in range(N):
    x = np.cos(t + i)
    y = np.sin(t + i)
    z = t / (2 + i*0.3)
    trajectories.append((x, y, z))

## Création figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.set_zlim(0, 10)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

## Couleurs automatiques
colors = plt.cm.viridis(np.linspace(0, 1, N))

## Listes des objets graphiques
points = []
#lines = []

for i in range(N):
    point, = ax.plot([], [], [], 'o', color=colors[i])
#    line, = ax.plot([], [], [], '-', color=colors[i], alpha=0.5)
    points.append(point)
    #lines.append(line)

def update(frame):
    for i in range(N):
        x, y, z = trajectories[i]

        ## Position actuelle (⚠ toujours mettre dans des listes)
        points[i].set_data([x[frame]], [y[frame]])
        points[i].set_3d_properties([z[frame]])

        ## Trajectoire passée
        #lines[i].set_data(x[:frame], y[:frame])
        #lines[i].set_3d_properties(z[:frame])

    return points #+ lines

ani = FuncAnimation(fig, update, frames=len(t), interval=40)

plt.show()