#The circumbinary Planets with a little Star Modification
from math import cos, sin, pi
import scipy as sci#Import matplotlib and associated modules for 3D and animations
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation as an
from IPython.display import HTML
from numpy import arange
from numpy import meshgrid
import numpy as np
import math as m
#Define universal gravitation constant
G=6.67408e-11 #N-m2/kg2#Reference quantities
m_nd=1.989e+30 #kg #mass of the sun
r_nd=5e+12 #m #distance between stars in Alpha Centauri
v_nd=30000 #m/s #relative velocity of earth around the sun
t_nd=79.91*365*24*3600*0.51 #s #orbital period of Alpha Centauri

#Net constants
K1=G*t_nd*m_nd/(r_nd**2*v_nd)
K2=v_nd*t_nd/r_nd
e=0.2
ep = (1-e**2)**(1/2)
#Define masses
m1=0.5 #Star A
m2=1 #Star B
m3 =0.0003 #Planet Dasvidania
a1 = m2/(m1+m2)
a2 = m1/(m1+m2)
#Define initial position vectors
r1=[a1+e*a1,0] #m
r2=[-a2-a2*e,0] #m
r3 = [2.25,0] #m
#Convert pos vectors to arrays
r1=sci.array(r1,dtype="float64")
r2=sci.array(r2,dtype="float64")
r3=sci.array(r3, dtype="float64")
#Find Centre of Mass

r_com=(m1*r1+m2*r2+m3*r3)/(m1+m2+m3)

#Define initial velocities
v3=[0,0.14] #m/s
 
#Convert velocity vectors to arrays
v3=sci.array(v3,dtype="float64")

#A function defining the equations of motion 
def ThreeBodyEquations(w,t,G,m1,m2,m3):
  
    x1 = -2*pi*(a1*sin(2*pi*t))
    y1 = -2*pi*(ep*a1*cos(2*pi*t))
    x2 = 2*pi*(a2*sin(2*pi*t))
    y2 = 2*pi*(ep*a2*cos(2*pi*t))
    r1=w[:2]
    r2=w[2:4]
    r3=w[4:6]
    v3=w[6:8]
    
    #Calculate magnitude or norm of vector  

    r12=sci.linalg.norm(r2-r1)
    r13=sci.linalg.norm(r3-r1)
    r23=sci.linalg.norm(r3-r2)
    
    #Differential Equations for planets
    
    dv3dt=K1*m1*(r1-r3)/r13**3 + K1*m2*(r2-r3)/r23**3
    dr1dt = [x1,y1]
    dr2dt = [x2,y2]
    dr3dt=K2*v3
    
    #Concatenation
    r12=sci.concatenate((dr1dt,dr2dt))
    r_derivs=sci.concatenate((r12,dr3dt))
    
    derivs=sci.concatenate((r_derivs,dv3dt))
    return derivs

 
#Package initial parameters
init_params=sci.array([r1,r2,r3,v3]) #create array of initial params
init_params=init_params.flatten() #flatten array to make it 1D
time_span=sci.linspace(0,25,500) #8 orbital periods and 500 points#Run the ODE solver
import scipy.integrate
#Solution to the system of differential equations 
three_body_sol=sci.integrate.odeint(ThreeBodyEquations,init_params,time_span,args=(G,m1,m2,m3))
r1_sol=three_body_sol[:,:2]
r2_sol=three_body_sol[:,2:4]
r3_sol=three_body_sol[:,4:6]
#The Center of mass: 
r1com_sol = r1_sol
r2com_sol = r2_sol
r3com_sol = r3_sol 

n=500
def x1(i):
      x1 = sum(r1com_sol[i:-n+i+1,0])
      return x1    
def y1(i):
      y1 = sum(r1com_sol[i:-n+i+1,1])
      return y1
def x2(i):
      x2 = sum(r2com_sol[i:-n+i+1,0])
      return x2    
def y2(i):
      y2 = sum(r2com_sol[i:-n+i+1,1])
      return y2  
delta = 0.025
xrange = arange(-15.0, 20.0, delta)
yrange = arange(-15.0, 20.0, delta)
X, Y = meshgrid(xrange,yrange)
def E(i):
    E = 1/(((X-x1(i))**2+(Y-y1(i))**2))+(1/16)/(pi*((X-x2(i))**2+(Y-y2(i))**2))
    return E
    
def update_points(i):
    # Updating orbits
    planet1_orbit.set_data(r1com_sol[:i, 0], r1com_sol[:i, 1])
    planet2_orbit.set_data(r2com_sol[:i, 0], r2com_sol[:i, 1]) 
    planet3_orbit.set_data(r3com_sol[:i, 0], r3com_sol[:i, 1])
 
    # Updating final points (The big circles)
    planet1_final.set_data(r1com_sol[i, 0], r1com_sol[i, 1])
    planet2_final.set_data(r2com_sol[i, 0], r2com_sol[i, 1])
    planet3_final.set_data(r3com_sol[i, 0], r3com_sol[i, 1])
 
#Create figure
fig, ax = plt.subplots()

#Plot the orbits
planet1_orbit, = ax.plot(r1com_sol[:,0], r1com_sol[:,1],   color="black", linewidth=0.13)
planet2_orbit, = ax.plot(r2com_sol[:,0], r2com_sol[:,1],  color="black", linewidth=0.13)
planet3_orbit, = ax.plot(r3com_sol[:,0], r3com_sol[:,1],  color="black", linewidth=0.4)

#Plot the final positions of the stars
# Plot final positions of the planets
planet1_final, = ax.plot(r1com_sol[0:1,0], r1com_sol[0:1,1],  color="darkblue", marker="o", markersize=6, label="Star A")
planet2_final, = ax.plot(r2com_sol[0:1,0], r2com_sol[0:1,1], color="gold", marker="o", markersize=6, label="Star B")
planet3_final, = ax.plot(r3com_sol[0:1,0], r3com_sol[0:1,1],  color="green", marker="o", markersize=5, label="Planet D")
for i in range(1, 500):
      plt.contourf(X, Y, E(i), [ 2*0.54/(4*m.pi), 2*1.1/(4*m.pi)],fill = True, colors= "green", alpha = 0.01)
 
ax.set_xlabel("$x$-coordinate",fontsize=14)
ax.set_ylabel("$y$-coordinate",fontsize=14)
plt.xlim(-6,6)
plt.ylim(-6,6)
plt.rcParams['figure.figsize'] = 12, 12
ax.set_title("Visualization of orbits of stars in a two-body system\n",fontsize=14)
ax.legend(loc="upper left",fontsize=14)
plt.grid()
anim = an.FuncAnimation(fig , update_points, frames = 500, interval = 50)
#writer = an.FFMpegWriter(fps=30,bitrate=1800)
#anim.save('ThreeBodyProblem.mp4', w
HTML(anim.to_html5_video())
