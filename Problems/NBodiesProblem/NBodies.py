import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import matplotlib.patches as patches
import numpy as np
import math
import random
import scipy.optimize
import matplotlib.animation as animation
    
def f(x, y, m, rs, masses):
    resx = 0
    resy = 0
    
    for i in range(0, len(rs), 4):
        x0 = rs[i]
        y0 = rs[i + 1]
       
        r = math.sqrt((x - x0) ** 2 + (y - y0) ** 2)
        
        if (r < 0.0001):
            continue
        
        fx = m * masses[int(i / 4)] * (x0 - x) / (r ** 3)
        fy = m * masses[int(i / 4)] * (y0 - y) / (r ** 3)
        
        FACTOR = 1000000000
        
        fx = min(fx, FACTOR)
        fy = min(fy, FACTOR)
        
        resx += fx
        resy += fy
        
    return resx, resy

def F(rs, masses):
    res = [0 for _ in range(len(rs))]
    
    for i in range(0, len(rs), 4):
        x = rs[i]
        y = rs[i + 1]
        vx = rs[i + 2]
        vy = rs[i + 3]
        
        m = masses[int(i / 4)]
        
        fx, fy = f(x, y, m, rs, masses)
        
        #
        # crazy values control
        #
        R = 300
        r2 = x**2 + y**2
        if r2 > R**2:
            rs[i + 0] = min(max(-R, rs[i + 0]), R)
            rs[i + 1] = min(max(-R, rs[i + 1]), R)
            rs[i + 2] = 0
            rs[i + 3] = 0
        #
        #
        #
        
        res[i] = vx
        res[i + 1] = vy
        res[i + 2] = fx / m
        res[i + 3] = fy / m

    return np.array(res)

def integrate(rs0, masses, t0, t1, n):
    h = (t1 - t0) / n
    t = np.arange(t0, t1, h)
    y = [np.array(range(len(rs0))) for _ in t]
    
    y[0] = np.array(rs0)
    
    for i in range(1, len(t)):
        k1 = F(y[i-1], masses)
        k2 = F(y[i-1] + 0.5*h*k1, masses)
        k3 = F(y[i-1] + 0.5*h*k2, masses)
        k4 = F(y[i-1] + h * k3, masses)

        y[i] = y[i-1] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
        

    return y


# default
drawEnergy = True
masses = [10, 1]
inits = [0, 0, 0, 0.1,      10, 0, 0, 1]
t0 = 0
t1 = 70
n = 1000

frameSpeed = 1
waitTime = 0
alp = 1


# solar system order
"""
masses = []#[100, 1, 1]
inits = []#[0, 0, 0, 0.1,      10, 0, 0, 1,       -10, 0, 0, -1]

N = 10

for i in range(N):
    masses.append(1)
    tetta = 2*math.pi * i / N
    r = 10
    x = r * math.cos(tetta)
    y = r * math.sin(tetta)
    v = 1#random.uniform(1, 2)
    dtetta = 0#random.uniform(0.1, 0.1)
    vx = -v * math.sin(tetta+dtetta)
    vy = v * math.cos(tetta+dtetta)
    
    inits.append(x)
    inits.append(y)
    inits.append(vx)
    inits.append(vy)

masses.append(10)
inits.append(0)
inits.append(0)
inits.append(0)
inits.append(0)

t0 = 0
t1 = 70
n = 1000

frameSpeed = 0.3
waitTime = 0#0.1
"""

# solar system chaos
"""
masses = []#[100, 1, 1]
inits = []#[0, 0, 0, 0.1,      10, 0, 0, 1,       -10, 0, 0, -1]

N = 10

for i in range(N):
    masses.append(1)
    tetta = 2*math.pi * i / N
    r = i * 10 + 70
    x = r * math.cos(tetta)
    y = r * math.sin(tetta)
    v = math.sqrt(100 / r)
    vx = -v * math.sin(tetta)
    vy = v * math.cos(tetta)
    
    inits.append(x)
    inits.append(y)
    inits.append(vx)
    inits.append(vy)

masses.append(100)
inits.append(0)
inits.append(0)
inits.append(0)
inits.append(0)

t0 = 0
t1 = 5000
n = 1000

frameSpeed = 3
waitTime = 0
alp = 0.3
"""

# calculations

res = integrate(inits, masses, t0, t1, n)


fig, axs, ax = None, None, None

if drawEnergy:
    fig, axs = plt.subplots(2)
    ax = axs[0]
else:
    fig, axs = plt.subplots(1)
    ax = axs

N = int(len(res[0]) / 4)

color = iter(cm.rainbow(np.linspace(0, 1, N)))

lines = []
points = []
ln = None

for i in range(N):
    xs = [res[j][i*4] for j in range(len(res))]
    ys = [res[j][i*4 + 1] for j in range(len(res))]
    
    c = next(color)
    sc = ax.scatter(xs[0], ys[0], color = c, s = 10 * masses[i])
    l, = ax.plot(xs, ys, label = "Body " + str(i + 1), color = c, alpha = alp)
    lines.append(l)
    points.append(sc)
        
def animate(i):
    i = min(int(i), len(res) - 1)
    if ln != None:
        ln.set_xdata([t[k] for k in range(min(int(i), len(t) - 1))])
        ln.set_ydata([vs2[k] for k in range(min(int(i), len(t) - 1))])
    for j in range(N):
        lines[j].set_xdata([res[k][j*4] for k in range(i)])
        lines[j].set_ydata([res[k][j*4 + 1] for k in range(i)])
        points[j].set_offsets([res[i][j*4], res[i][j*4 + 1]])
    return [lines, points]

#ani = animation.FuncAnimation(fig, animate, interval=1, frames = np.arange(0, (1+waitTime)*len(res), (1+waitTime)*len(res) * frameSpeed / (t1-t0)))

ax.legend()
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.grid()
ax.legend()

# draw energy
if drawEnergy:
    vxs = [res[j][0*4 + 2] for j in range(len(res))]
    vys = [res[j][0*4 + 3] for j in range(len(res))]
    t = np.arange(t0, t1, (t1 - t0) / n)
    vs2 = [vxs[i]**2 + vys[i]**2 for i in range(len(vxs))]

    ln, = axs[1].plot(t, vs2, label = "Body 1 kinetic energy")
    
    def animate2(i):
        ln.set_xdata([t[k] for k in range(min(int(i), len(t) - 1))])
        ln.set_ydata([vs2[k] for k in range(min(int(i), len(t) - 1))])
        return ln

    #aniE = animation.FuncAnimation(fig, animate2, interval=1, frames = np.arange(0, (1+waitTime)*len(res), (1+waitTime)*len(res) * frameSpeed / (t1-t0)))
    
    axs[1].legend()
    axs[1].set_xlabel("t")
    axs[1].set_ylabel("E")
    axs[1].grid()
    axs[1].legend()
    
ani = animation.FuncAnimation(fig, animate, interval=1, frames = np.arange(0, (1+waitTime)*len(res), (1+waitTime)*len(res) * frameSpeed / (t1-t0)))
    
plt.show()
        
ani.save('animation.gif', fps=25)
    


