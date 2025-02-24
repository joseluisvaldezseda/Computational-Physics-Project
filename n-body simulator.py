import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation as an
from IPython.display import HTML
#n-body equations.
n = int(input("How many stars do you want the most?: "))
def n_body_equations(x, t, n):
    
   # Unpacking positions
   r = [] 
   for i in range(0,n):
        r.append(x[3*i:3*(i+1)])
    
   # Unpacking velocities
   v = []
   for i in range(n, 2*n):
        v.append(x[3*i:3*(i+1)])
 
    # Constants
   G = 6.67408e-10
    
   out = np.zeros(3*2*n)
     
   # out         = np.zeros(18)
   for i in range(0,n):
        out[3*i:3*(i+1)] = K2*v[i]
    
   k = 1
   sums = []
   for i in range(n):
        s=0
        j = 0
        while j<n:
            if j != i:
                term = K1*(r[j]-r[i])/(np.linalg.norm(r[j]-r[i]))
                s = s + term
            j = j+1 
        sums.append(s)
    
   for i in range(n,2*n):
        out[3*i:3*(i+1)]= sums[i-n]
    
   return out
 # Reference quantities
m_nd = 1.989e+30                # mass of sun                               [kg]
r_nd = 5.326e+12                # Distance between stars in Alpha Centauri  [m]
v_nd = 30000                    # Relative velocity of earth around the sun [m/s]
t_nd = 79.91*365*24*3600*0.51   # Orbital period of Alpha Centauri          [s]
 
# Constants
G = 6.67408e-10
K1 = G*t_nd*m_nd/((r_nd ** 2)*v_nd)
K2 = v_nd*t_nd/r_nd
 
## Parameters and initial conditions
# Planet masses: x times the mass of the sun   
    
import random
# Define initial positions
r0 = []
for i in range(n):
    r0.append(np.array([random.uniform(-40,40) ,random.uniform(-40,40),0.0]))
Tr = tuple(r0)
 
# Define initial velocities
v0= []
for i in range(n):
    v0.append(np.array([random.uniform(-0.5,0.5),random.uniform(-0.5,0.5),random.uniform(-0.1,0.1)]))
Tv = tuple(v0)
 
# Initial conditions
#x0          = np.zeros(18)
x0 = np.zeros(3*2*n)
 
x0[0:3*n]         = np.concatenate(Tr)
 
x0[3*n: 3*(2*n)]  = np.concatenate(Tv)
 
## Integrating ODE
tf = 5
c = 500
t = np.linspace(0, tf, num = c)
n_body_sol = odeint(n_body_equations, x0, t, args=(n,))
r_sol = []
for i in range(n):
    r_sol.append(n_body_sol[:, 3*i:3*(i+1)])
 
## Plotting
# Animation function
 
def update_points(i):
    # Updating orbits
    j = 0
    for r in r_sol:
        j = j+1
        globals()["planet"+str(j)+"_orbit"].set_data(r[:i, 0], r[:i, 1])
        globals()["planet"+str(j)+"_orbit"].set_3d_properties(r[:i, 2], 'z')
  
    # Updating points  
    k = 0
    for r in r_sol:
        k = k+1
        globals()["planet"+str(k)+"_final"].set_data(r[i, 0], r[i, 1])
        globals()["planet"+str(k)+"_final"].set_3d_properties(r[i, 2], 'z')
 
# Create figure
fig = plt.figure(figsize=(8, 8))
 
# Create 3D axes
ax = fig.add_subplot(111, projection="3d")
 
# Initial orbits
i = 0
for r in r_sol:
    i = i+1
    globals()["planet"+str(i)+"_orbit"], = ax.plot(r[:,0], r[:,1], r[:,2], color="darkblue", linewidth=0)
j = 0
for r in r_sol:
    j = j+1
    globals()["planet"+str(j)+"_final"], = ax.plot(r[0:1,0], r[0:1,1], r[0:1,2], color="darkblue", marker="o", markersize=3)

# Misc. plot details
ax.set_xlabel("x-coordinate", fontsize=14)
ax.set_ylabel("y-coordinate", fontsize=14)
#ax.set_zlabel("z-coordinate", fontsize=7)
ax.set_title("n body animation", fontsize=14)
ax.legend(loc="upper left", fontsize=8)
ax.set_xlim(-45,45)
ax.set_ylim(-45,45)
ax.set_zlim(-2,2)
# Create animation
 
# Export animation
anim = an.FuncAnimation(fig , update_points, frames = 500, interval = 50)
#writer = an.FFMpegWriter(fps=30,bitrate=1800)
#anim.save('ThreeBodyProblem.mp4', writer=writer)
HTML(anim.to_html5_video())
