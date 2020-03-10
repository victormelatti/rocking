#theta in function of time for accelerogram ground motion
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import interpolate

b=0.38 #mezza base
h=1.6 #mezza altezza
zeta=h/b
vol=2*b*2*h
density=203
g=9.81
m=vol*density
R=np.sqrt((b)**2+(h)**2)
I=(4/3)*m*R**2
p_quadro=m*g*R/I
p=np.sqrt(p_quadro)

#u_teta
B=2*b
l_h=1
k_n=6e6+6e6*0.5 #Pa
f_m=32739+32739*0.5 #Pa 
teta_j0=2*m*g/(k_n*B**2*l_h)
a_c=2*m*g/(f_m*l_h)
teta_jc=f_m/(k_n*a_c)


list1=[]
list2=[]
def u_teta_dyn(theta,list10,list20):
    global list1
    global list2
    u_theta=0
    if theta<teta_j0 and len(list1)==0:
        u_theta=(B/2-(B**3*k_n*l_h*theta/(12*m*g)))
    if theta <=teta_jc and theta >teta_j0 and len(list2)==0 or (len(list1)!=0 and len(list2)==0): 
        u_theta=(1/3)*np.sqrt(2*m*g/(k_n*l_h*theta))
        list1.append(theta)
    if theta>teta_jc or (len(list1)!=0 and len(list2)!=0):
        u_theta=0.5*(m*g/(f_m*l_h)+(f_m**3*l_h)/(12*m*g*k_n**2*theta**2))
        list2.append(theta)
    R_theta=np.sqrt((R*np.cos(alpha))**2+(R*np.sin(alpha)-u_theta)**2)
    I_theta=m/3*R**2+m*R_theta**2
    p_theta=np.sqrt(m*g*R/I_theta)
    return u_theta, p_theta

def f(y, t, params): 
    theta,omega = y
    p_quadro, alpha, interpol,g, u_teta_dyn, list1, list2 = params
    derivs = [omega,
              -u_teta_dyn(theta,list1,list2)[1]**2*(np.sin(alpha-theta)-u_teta_dyn(theta,list1,list2)[0]/R)-p_quadro*interpol(t)/g*np.cos(alpha-theta)]
    return derivs
# Make time array for solution
tStop = 8.004
tInc = 0.005
t = np.arange(0., tStop, tInc)
# Parameters
alpha=np.arctan(b/h)

import xlrd
wb=xlrd.open_workbook('D:\\PHD\\TEACHING\\ANALISI DINAMICA NON LINEARE\\acceleration_laquila.xlsx')
s=wb.sheet_by_name('acceleration strong')
s.nrows
time=np.arange(0., tStop, tInc)
#time=time[6000:10000]
acc=[s.cell_value(i,4)*0.01 for i in range(s.nrows)]
#acc=acc[6000:10000]
interpol=interpolate.interp1d(time,acc,fill_value="extrapolate")
#plt.plot(time,acc)
#plt.plot(t,interpol(t))

fig, ax = plt.subplots()
axes = plt.gca()
ax.plot(t,interpol(t))
ax.set(xlabel='Time [m]',ylabel='rotations [degree]')
#axes.set_ylim([-10,10]) 
axes.set_xlim([0,10]) 

# Initial values
theta0 = 5*np.pi/180     # initial angular displacement RADIANTI
omega0 = 0.0     # initial angular velocity
params = [p_quadro, alpha, interpol,g, u_teta_dyn, list1, list2]
y0 = [theta0, omega0]

# Call the ODE solver
psoln = odeint(f, y0, t, args=(params,))

# Plot theta as a function of time

fig, ax = plt.subplots()
axes = plt.gca()
ax.plot(t, psoln[:,0]*180/np.pi)
ax.set(xlabel='Time [m]',ylabel='rotations [degree]')
#axes.set_ylim([-10,10]) 
#axes.set_xlim([0,2])   
#%%
#AGGIUNGO L'URTO

#e_1s=(1-1.5*np.sin(alpha)**2)**2*(1-1.5*np.cos(alpha)**2)
e_1s=1.05*(1-2*m*R**2/I*np.sin(alpha)**2)**2*(1-2*m*R**2/I*np.cos(alpha)**2)
rotations = []
velocities = []
ts = [t]
numero_di_volte=7
for ndv in range(numero_di_volte):
    rotations.append([])
    velocities.append([])
    for i in range(len(psoln[:,0])):
        if psoln[i,0]>=0:
            rotations[-1]=rotations[-1]+[psoln[i,0]*180/np.pi]
            velocities[-1]=velocities[-1]+[psoln[i,1]]
        else:
            break
    
    fig, ax = plt.subplots()
    axes = plt.gca()
    ax.grid()
    ax.plot(ts[-1][:len(rotations[-1])],rotations[-1],'r-',linewidth=1)
    ax.set(xlabel='time (sec)',ylabel='$\Theta$ (deg)',title='singular dynamic of $\Theta$')
    #primo urto
    r=e_1s
    velocity_urto_meno=velocities[-1][-1]
    velocity_urto_piu=r*velocity_urto_meno

    ts.append([])
    ts[-1] = np.arange(ts[-2][len(rotations[-1])], tStop, tInc)
    y_curr=[0, velocity_urto_piu]
    interpol=interpolate.interp1d(time[len(rotations[-1]):],acc[len(rotations[-1]):],fill_value="extrapolate")
    #interpol2=interpolate.interp1d(time[len(rotations[-1]):],acc[len(rotations[-1]):])

    params = [p_quadro, alpha, interpol,g, u_teta_dyn, list1, list2]
    
    psoln = odeint(f, y_curr, ts[-1], args=(params,))
#np.array(rotations[0]).min()    
#METTO INSIEME
t_total = []
rotation_total = []
for ndv in range(numero_di_volte):
    t_total += list(ts[ndv][:len(rotations[ndv])])
    rotation_total += rotations[ndv]

fig, ax = plt.subplots()
axes = plt.gca()
ax.grid()
ax.plot(t_total,rotation_total,'r-',linewidth=3)
ax.plot(time,acc,linewidth=0.5)
ax.set(xlabel='time (sec)',ylabel='$\Theta$ (deg)',title='complete dynamic of $\Theta$')
axes = plt.gca()
axes.set_xlim([0,10])
#axes.set_ylim([-1,10])

#fig.savefig("Plot theta as a function of time_WITH IMPACT DISSIPATION_igor_0.5.png", dpi = 500, bbox_inches='tight')
#%%




