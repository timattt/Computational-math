import numpy as np
import matplotlib.pyplot as plt


def lorenz_system(x, y, z, r, b=10, s=6):
    x_dot = b * (y - x)
    y_dot = r * x - y - x * z
    z_dot = x * y - s * z
    return x_dot, y_dot, z_dot


dr = 0.1
r = np.arange(40, 200, dr)
dt = 0.001
t = np.arange(0, 10, dt)

xs = np.empty(len(t) + 1)
ys = np.empty(len(t) + 1)
zs = np.empty(len(t) + 1)

xs[0], ys[0], zs[0] = (1, 1, 1)

r_maxes = []
z_maxes = []
r_mins = []
z_mins = []

for R in r:
    for i in range(len(t)):
        x_dot, y_dot, z_dot = lorenz_system(xs[i], ys[i], zs[i], R)
        xs[i + 1] = xs[i] + (x_dot * dt)
        ys[i + 1] = ys[i] + (y_dot * dt)
        zs[i + 1] = zs[i] + (z_dot * dt)
    for i in range(1, len(zs) - 1):
        if zs[i - 1] < zs[i] and zs[i] > zs[i + 1]:
            r_maxes.append(R)
            z_maxes.append(zs[i])
        elif zs[i - 1] > zs[i] and zs[i] < zs[i + 1]:
            r_mins.append(R)
            z_mins.append(zs[i])

    # "use final values from one run as initial conditions for the next to stay near the attractor"
    xs[0], ys[0], zs[0] = xs[i], ys[i], zs[i]

fig, ax = plt.subplots()

ax.scatter(r_maxes, z_maxes, color="black", s=0.5, alpha=0.2, label = "maximums")
ax.scatter(r_mins, z_mins, color="red", s=0.5, alpha=0.2, label = "minimums")

plt.xlim(0, 200)
plt.ylim(0, 400)

ax.grid()
ax.set_xlabel("r parameter")
ax.set_ylabel("z value")

plt.show()