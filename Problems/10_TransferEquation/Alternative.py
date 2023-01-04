import numpy as np
from matplotlib import animation
from scipy import interpolate
from numpy import where
from math import sin
import math

import matplotlib; matplotlib.use('Qt4Agg')
import matplotlib.pylab as plt
plt.get_current_fig_manager().window.raise_()


LNWDT=2; FNT=15
plt.rcParams['lines.linewidth'] = LNWDT; plt.rcParams['font.size'] = FNT


init_func=4   # Select stair case function (0) or sin^2 function (1)

# function defining the initial condition
if (init_func==0):
    def f(x):
        f = np.zeros_like(x)
        f[np.where(x <= 0.1)] = 1.0
        return f
elif(init_func==1):
    def f(x):
        f = np.zeros_like(x)
        x_left = 0.25
        x_right = 0.75
        xm = (x_right+x_left)/2.0
        f = where((x>x_left) & (x<x_right), np.exp(-100*(x-xm)**2),f) 
        return f
elif(init_func==2):
    def f(x):
        f = np.zeros_like(x)
        f[np.where((x > 0.1) )] = 1.0
        return f
elif(init_func==3):
    def f(x):
        f = np.zeros_like(x)
        f[np.where((x > 0.3) &(x < 0.7) )] = 1.0
        return f
elif(init_func==4):
    def f(x):
        f = np.zeros_like(x)
        f = where((x>0) & (x<1), np.exp(x),f) 
        return f

def ftbs(u): # forward time backward space
    u[1:-1] = (1-c)*u[1:-1] + c*u[:-2]
    return u[1:-1]

def implicit(u): # forward time backward space
    u[1:-1] = (u[:-2] +u[2:])/2.0 -  c*(u[2:] - u[:-2])/2.0
    return u[1:-1]

# Lax-Wendroff
def lax_wendroff(u): 
    u[1:-1] = c/2.0*(1+c)*u[:-2] + (1-c**2)*u[1:-1] - c/2.0*(1-c)*u[2:]
    return u[1:-1]

def experimental(u):
    cp = u.copy()
    cp[1:-1] = (c**2 - c) / 2 * cp[:-2] + cp[1:-1]*(c-c**2)+cp[2:]*(c-c**2)/2
    
    u[1:-1] = (1-c)*u[1:-1] + c*u[:-2]
    
    i1 = 0
    i2 = 0
    
    for i in range(len(cp)-1):
        if (u[i] == 0 and u[i+1] != 0) or (u[i] != 0 and u[i+1] == 0):
            if i1 == 0:
                i1 = i
            else:
                i2 = i
                
    r = 1
    for i in range(len(cp)-1):
        if math.fabs(i - i1) < r or math.fabs(i - i2) < r:
            u[i] += cp[i]
            
    return u[1:-1]

# Constants and parameters
a = 1.0 # wave speed
tmin, tmax = 0.0, 1.0 # start and stop time of simulation
xmin, xmax = 0.0, 2.0 # start and end of spatial domain
Nx = 80 # number of spatial points
c = 0.9 # courant number, need c<=1 for stability


# Discretize
x = np.linspace(xmin, xmax, Nx+1) # discretization of space
dx = float((xmax-xmin)/Nx) # spatial step size
dt = c/a*dx # stable time step calculated from stability requirement
Nt = int((tmax-tmin)/dt) # number of time steps
time = np.linspace(tmin, tmax, Nt) # discretization of time

# solve from tmin to tmax

solvers = [ftbs, implicit, lax_wendroff]
#solvers = [ftbs,lax_wendroff,macCormack]
#solvers = [ftbs,lax_wendroff]
solvers = [experimental]

u_solutions=np.zeros((len(solvers),len(time),len(x)))
uanalytical = np.zeros((len(time), len(x))) # holds the analytical solution


    
for k, solver in enumerate(solvers): # Solve for all solvers in list
    u = f(x)
    un = np.zeros((len(time), len(x))) # holds the numerical solution

    for i, t in enumerate(time[1:]):
        
        if k==0:
            uanalytical[i,:] = f(x-a*t) # compute analytical solution for this time step
            
        u_bc = interpolate.interp1d(x[-2:], u[-2:]) # interplate at right bndry
        
        u[1:-1] = solver(u[:]) # calculate numerical solution of interior
        u[-1] = u_bc(x[-1] - a*dt) # interpolate along a characteristic to find the boundary value
        
        un[i,:] = u[:] # storing the solution for plotting
    
    u_solutions[k,:,:] = un



### Animation 
 
# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(xmin,xmax), ylim=(np.min(un), np.max(un)*1.1))

lines=[]     # list for plot lines for solvers and analytical solutions
legends=[]   # list for legends for solvers and analytical solutions

for solver in solvers:
    line, = ax.plot([], [])
    lines.append(line)
    legends.append(solver.__name__)

line, = ax.plot([], []) #add extra plot line for analytical solution
lines.append(line)
legends.append('Analytical')

plt.xlabel('x-coordinate [-]')
plt.ylabel('Amplitude [-]')
plt.legend(legends, loc=3, frameon=False)
 
# initialization function: plot the background of each frame
def init():
    for line in lines:
        line.set_data([], [])
    return lines,

# animation function.  This is called sequentially
def animate(i):
    for k, line in enumerate(lines):
        if (k==0):
            line.set_data(x, un[i,:])
        else:
            line.set_data(x, uanalytical[i,:])
    return lines,

def animate_alt(i):
    for k, line in enumerate(lines):
        if (k==len(lines)-1):
            line.set_data(x, uanalytical[i,:])
        else:
            line.set_data(x, u_solutions[k,i,:])
    return lines,

 
# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate_alt, init_func=init, frames=Nt, interval=100, blit=False)
anim.save('animation.gif', fps=25)
 
plt.show()