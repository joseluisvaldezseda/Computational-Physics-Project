#Import scipy
import scipy as sci#Import matplotlib and associated modules for 3D and animations
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

#Define universal gravitation constant
G=6.67408e-11 #N-m2/kg2#Reference quantities
m_nd=1.989e+30 #kg #mass of the sun
r_nd=5.326e+12 #m #distance between stars in Alpha Centauri
v_nd=30000 #m/s #relative velocity of earth around the sun
t_nd=79.91*365*24*3600*0.51 #s #orbital period of Alpha Centauri#Net constants
K1=G*t_nd*m_nd/(r_nd**2*v_nd)
K2=v_nd*t_nd/r_nd

#Define masses
m1=1.1 #Alpha Centauri A
m2=0.0004 #Alpha Centauri B#Define initial position vectors
r1=[-0.7,0,0] #m
r2=[0.7,0,0] #m#Convert pos vectors to arrays
r1=sci.array(r1,dtype="float64")
r2=sci.array(r2,dtype="float64")#Find Centre of Mass
r_com=(m1*r1+m2*r2)/(m1+m2)#Define initial velocities
v1=[0,0,0] #m/s
v2=[0,0.1,0] #m/s#Convert velocity vectors to arrays
v1=sci.array(v1,dtype="float64")
v2=sci.array(v2,dtype="float64")#Find velocity of COM
v_com=(m1*v1+m2*v2)/(m1+m2)

#A function defining the equations of motion 
def TwoBodyEquations(w,t,G,m1,m2):
    r1=w[:3]
    r2=w[3:6]
    v1=w[6:9]
    v2=w[9:12]    
    r=sci.linalg.norm(r2-r1)
    dv1bydt=K1*m2*(r2-r1)/r**3     #Calculate magnitude or norm of vector    dv1bydt=K1*m2*(r2-r1)/r**3
    dv2bydt=K1*m1*(r1-r2)/r**3
    dr1bydt=K2*v1
    dr2bydt=K2*v2    
    r_derivs=sci.concatenate((dr1bydt,dr2bydt))
    derivs=sci.concatenate((r_derivs,dv1bydt,dv2bydt))
    return derivs

#Package initial parameters
#Package initial parameters
init_params=sci.array([r1,r2,v1,v2]) #create array of initial params
init_params=init_params.flatten() #flatten array to make it 1D
time_span=sci.linspace(0,20, 600) #8 orbital periods and 500 points#Run the ODE solver
import scipy.integrate

two_body_sol=sci.integrate.odeint(TwoBodyEquations,init_params,time_span,args=(G,m1,m2))
r1e_sol=two_body_sol[:,:3]
r2e_sol=two_body_sol[:,3:6]
r1_sol = r1e_sol -r_com
r2_sol = r2e_sol -r_com
#Create figure
fig=plt.figure(figsize=(15,15))#Create 3D axes
ax=fig.add_subplot(111)#Plot the orbits
plt.xlim(-6,6)
plt.ylim(-6,6)
ax.plot(r1_sol[:,0],r1_sol[:,1] ,color="darkblue")
ax.plot(r2_sol[:,0],r2_sol[:,1] ,color="black",lw=0.3)#Plot the final positions of the stars
ax.scatter(r1_sol[-1,0],r1_sol[-1,1], color="gold",marker="o",s=100,label="Alpha Centauri A")
ax.scatter(r2_sol[-1,0],r2_sol[-1,1], color="darkblue",marker="o",s=50,label="Alpha Centauri B")#Add a few more bells and whistles
ax.set_xlabel("x-coordinate",fontsize=14)
ax.set_ylabel("y-coordinate",fontsize=14)
plt.grid()
#How do I make an animation of this?
