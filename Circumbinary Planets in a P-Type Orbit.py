#Import scipy
import scipy as sci#Import matplotlib and associated modules for 3D and animations
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation as an
from IPython.display import HTML
#Define universal gravitation constant
G=6.67408e-11 #N-m2/kg2#Reference quantities
m_nd=1.989e+30 #kg #mass of the sun
r_nd=5.326e+12 #m #distance between stars in Alpha Centauri
v_nd=30000 #m/s #relative velocity of earth around the sun
t_nd=79.91*365*24*3600*0.51 #s #orbital period of Alpha Centauri

#Net constants
K1=G*t_nd*m_nd/(r_nd**2*v_nd)
K2=v_nd*t_nd/r_nd

#Define masses
m1=0.5 #Star A
m2=1 #Star B
m3 =0.0003 #Planet C
#Define initial position vectors
r1=[-0.5,0,0] #m
r2=[0.5,0,0] #m
r3 = [2,0,0] #m
#Convert pos vectors to arrays
r1=sci.array(r1,dtype="float64")
r2=sci.array(r2,dtype="float64")
r3=sci.array(r3, dtype="float64")
#Find Centre of Mass

r_com=(m1*r1+m2*r2+m3*r3)/(m1+m2+m3)

#Define initial velocities
v1=[0.01,0.01,0] #m/s
v2=[-0.05,-0.1,0] #m/s
v3=[0,0.1,0]
#Convert velocity vectors to arrays
v1=sci.array(v1,dtype="float64")
v2=sci.array(v2,dtype="float64")
v3=sci.array(v3,dtype="float64")


#A function defining the equations of motion 
def ThreeBodyEquations(w,t,G,m1,m2,m3):
    r1=w[:3]
    r2=w[3:6]
    r3=w[6:9]
    v1=w[9:12]
    v2=w[12:15]    
    v3=w[15:18]
    
    #Calculate magnitude or norm of vector  

    r12=sci.linalg.norm(r2-r1)
    r13=sci.linalg.norm(r3-r1)
    r23=sci.linalg.norm(r3-r2)
    
    #Differential Equations for planets
    
    dv1dt=K1*m2*(r2-r1)/r12**3 + K1*m3*(r3-r1)/r13**3    
    dv2dt=K1*m1*(r1-r2)/r12**3 + K1*m3*(r3-r2)/r23**3
    dv3dt=K1*m1*(r1-r3)/r13**3 + K1*m2*(r2-r3)/r23**3
    
    dr1dt=K2*v1
    dr2dt=K2*v2
    dr3dt=K2*v3
    
    #Concatenation
    r12_derivs=sci.concatenate((dr1dt,dr2dt))
    r_derivs=sci.concatenate((r12_derivs,dr3dt))
    v12_derivs=sci.concatenate((dv1dt,dv2dt))
    v_derivs=sci.concatenate((v12_derivs,dv3dt))
    
    derivs=sci.concatenate((r_derivs,v_derivs))
    return derivs

 
#Package initial parameters
init_params=sci.array([r1,r2,r3,v1,v2,v3]) #create array of initial params
init_params=init_params.flatten() #flatten array to make it 1D
time_span=sci.linspace(0,20,500) #8 orbital periods and 500 points#Run the ODE solver
import scipy.integrate
#Solution to the system of differential equations 
three_body_sol=sci.integrate.odeint(ThreeBodyEquations,init_params,time_span,args=(G,m1,m2,m3))
r1_sol=three_body_sol[:,:3]
r2_sol=three_body_sol[:,3:6]
r3_sol=three_body_sol[:,6:9]
#The Center of mass:
rcom_sol = (m1*r1_sol + m2*r2_sol + m3*r3_sol)/(m1+m2)
r1com_sol = r1_sol-rcom_sol
r2com_sol = r2_sol-rcom_sol
r3com_sol = r3_sol-rcom_sol

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
planet1_final, = ax.plot(r1com_sol[0:1,0], r1com_sol[0:1,1],  color="darkblue", marker="o", markersize=8, label="Star A")
planet2_final, = ax.plot(r2com_sol[0:1,0], r2com_sol[0:1,1], color="gold", marker="o", markersize=10, label="Star B")
planet3_final, = ax.plot(r3com_sol[0:1,0], r3com_sol[0:1,1],  color="green", marker="o", markersize=5, label="Planet D")
 
ax.set_xlabel("$x$-coordinate",fontsize=14)
ax.set_ylabel("$y$-coordinate",fontsize=14)
plt.xlim(-4,4)
plt.ylim(-4,4)
plt.rcParams['figure.figsize'] = 10, 10
ax.set_title("Dasvidania\n",fontsize=14)
ax.legend(loc="upper left",fontsize=14)
plt.grid()
anim = an.FuncAnimation(fig , update_points, frames = 500, interval = 50)
#writer = an.FFMpegWriter(fps=30,bitrate=1800)
#anim.save('ThreeBodyProblem.mp4', writer=writer)
HTML(anim.to_html5_video())
 
