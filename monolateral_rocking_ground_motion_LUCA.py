#ground motion
#theta in function of time for accelerogram ground motion
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import xlrd
from scipy import interpolate
import copy

b=0.3 #mezza base
h=2.4 #mezza altezza
H=2*h
zeta=h/b
vol=2*b*2*h
density=203
g=9.81
m=vol*density
R=np.sqrt((b)**2+(h)**2)
I=(4/3)*m*R**2
p_quadro=m*g*R/I
p=np.sqrt(p_quadro)
alpha=np.arctan(b/h)

# Make time array for solution
tStop = 8.000
tInc = 0.005
t = np.linspace(0.0, tStop, int(tStop/tInc+1))


#seismic action
wb=xlrd.open_workbook('C:\\Users\\victor\\Desktop\\PHD\\rocking\\acceleration_laquila.xlsx')
s=wb.sheet_by_name('acceleration strong')
s.nrows
#time=np.linspace(0., tStop, int(tStop/tInc)+1)
#time=time[6000:10000]
acc=[s.cell_value(i,4)*0.01 for i in range(s.nrows)]

#extrapolate creates values after the last one
interpol=interpolate.interp1d(t,acc,fill_value="extrapolate")

def f(y, t, params): 
    theta,omega = y
    p_quadro, alpha, interpol,g = params
    
    derivs = [omega,
              -p_quadro*(np.sin(alpha-theta))+p_quadro*interpol(t)*np.cos(alpha-theta)/g]
    return derivs      
# Initial values
theta0 = 0*np.pi/180     # initial angular displacement RADIANTI
omega0 = 0.0  # initial angular velocity
params = [p_quadro, alpha, interpol,g]
y_curr = [theta0, omega0]

e_1s=1.05*(1-2*m*R**2/I*np.sin(alpha)**2)**2*(1-2*m*R**2/I*np.cos(alpha)**2)
rotations = []
velocities = []
ts = [t]
psoln = None
activation=0
#%%tests
numero_di_volte=90
activation_acceleration=[]
y_curr_list=[]


for ndv in range(numero_di_volte):    
    rotations.append([])
    velocities.append([])
    
    #velocity_urto_piu=0
    
    if psoln is not None:
        if psoln[0][1,0]>0:
            activation=1
            for i in range(len(psoln[0][:,0])):
                if psoln[0][i,0]>=0:
                    rotations[-1]=rotations[-1]+[psoln[0][i,0]*180/np.pi]
                    velocities[-1]=velocities[-1]+[psoln[0][i,1]]
                    r=e_1s 
                    velocity_urto_meno=velocities[-1][-1]
                    velocity_urto_piu=r*velocity_urto_meno
                    #activation_acceleration+=[acc[acc_new_start]]
                    counter=0
                else:
                    break        
        else:
            activation=0
            velocity_urto_piu=0
    
    if psoln is None:
        velocity_urto_piu=0

    y_curr=[0, velocity_urto_piu]
    y_curr_list+=[y_curr]
    if activation>0:
        ts.append([])
        new_start=ts[-2][len(rotations[-1])]
        ts[-1] = np.linspace(new_start, tStop, int(tStop/tInc+1)-int(new_start/tInc))
        #ts[-1] = np.arange(new_start, tStop+0.004, tInc)
        interpol=interpolate.interp1d(ts[-1],acc[int(new_start/tInc):],fill_value="extrapolate")
    if activation==0: # non ci sono stati valori buoni nel psoln 
        ts.append([])
        new_start=ts[-2][len(rotations[-1])]+ tInc
        ts[-1] = np.linspace(new_start, tStop, int(tStop/tInc+1)-int(new_start/tInc))
        #ts[-1] = np.arange(new_start + counter*tInc, tStop+0.004, tInc)
        interpol=interpolate.interp1d(ts[-1],acc[int(new_start/tInc):],fill_value="extrapolate")
    print(activation, new_start)
    params = [p_quadro, alpha, interpol,g]
    psoln = odeint(f, y_curr, ts[-1], args=(params,),full_output=True) 
#%%
#check lengths
    
len(ts[-1])
len(acc[int(new_start/tInc)+1:]) 
len(acc[int(new_start/tInc):]) 
#%%
        if new_start==0: # se il tratto precedente é stato senza valori buoni
            ts[-1] = np.arange(new_start + counter*tInc, tStop, tInc)
            interpol=interpolate.interp1d(t[new_start+counter:],acc[new_start+counter:],fill_value="extrapolate")
        else: # se il tratto precedente ha avuto valori positivi di rotazione ma poi nn ha fatto un balzo
            ts[-1] = np.arange(ts[-2][new_start + counter*tInc], tStop, tInc)
    print(activation)
    params = [p_quadro, alpha, interpol,g]
    psoln = odeint(f, y_curr, ts[-1], args=(params,),full_output=True)  

#%%
#np.array(rotations[0]).min()    
#METTO INSIEME
t_total = []
rotation_total = []
for ndv in range(numero_di_volte):
    t_total += list(ts[ndv][:len(rotations[ndv])])
    rotation_total += rotations[ndv]
initial_part_time=list(np.linspace(0.,t_total[0],int(t_total[0]/tInc+1)))
initial_part_rotations=[0]*len(initial_part_time)

t_total=initial_part_time+t_total
rotation_total=initial_part_rotations+rotation_total
interpol=interpolate.interp1d(t,acc)

fig, ax = plt.subplots()
axes = plt.gca()
ax.grid()
ax.plot(t_total,rotation_total,'r-',linewidth=3)
ax.plot(t_total,interpol(t_total))
ax.set(xlabel='time (sec)',ylabel='$\Theta$ (deg)',title='complete dynamic of $\Theta$')
axes = plt.gca()
axes.set_xlim([0,4])
#axes.set_ylim([-1,10])
#fig.savefig("Plot theta as a function of time_WITH IMPACT DISSIPATION_igor_0.5.png", dpi = 500, bbox_inches='tight')
#%%checking artificial ground motion: to test if a second acceleretion equal to the first one will
#produce the same acceleration
#theta in function of time for accelerogram ground motion
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import xlrd
from scipy import interpolate

b=0.3 #mezza base
h=2.4 #mezza altezza
H=2*h
zeta=h/b
vol=2*b*2*h
density=203
g=9.81
m=vol*density
R=np.sqrt((b)**2+(h)**2)
I=(4/3)*m*R**2
p_quadro=m*g*R/I
p=np.sqrt(p_quadro)
alpha=np.arctan(b/h)

# Make time array for solution
tStop = 3.000
tInc = 0.005
t = np.linspace(0., tStop, int(tStop/tInc)+1)

#seismic action
wb=xlrd.open_workbook('C:\\Users\\victor\\Desktop\\PHD\\rocking\\acceleration_laquila_checking.xlsx')
s=wb.sheet_by_name('acceleration strong')
s.nrows
#time=np.linspace(0., tStop, int(tStop/tInc)+1)
#time=time[6000:10000]
acc=[s.cell_value(i,4)*0.01 for i in range(s.nrows)]

#extrapolate creates values after the last one
interpol=interpolate.interp1d(t,acc,fill_value="extrapolate")

def f(y, t, params): 
    theta,omega = y
    p_quadro, alpha, interpol,g = params
    
    derivs = [omega,
              -p_quadro*(np.sin(alpha-theta))+p_quadro*interpol(t)*np.cos(alpha-theta)/g]
    return derivs      
# Initial values
theta0 = 0*np.pi/180     # initial angular displacement RADIANTI
omega0 = 0.0  # initial angular velocity
params = [p_quadro, alpha, interpol,g]
y_curr = [theta0, omega0]

e_1s=1.05*(1-2*m*R**2/I*np.sin(alpha)**2)**2*(1-2*m*R**2/I*np.cos(alpha)**2)
rotations = []
velocities = []
ts = [t]
psoln = None
#%%tests
numero_di_volte=10
psoln = None

for ndv in range(numero_di_volte): 
    rotations.append([])
    velocities.append([])
    activation=0
    
    if psoln is not None:
        for i in range(len(psoln[0][:,0])):
            if psoln[0][i,0]>=0:
                activation=1
                
                rotations[-1]=rotations[-1]+[psoln[0][i,0]*180/np.pi]
                velocities[-1]=velocities[-1]+[psoln[0][i,1]]
            
            else:
                break   
    
    r=e_1s if psoln is not None else 0
    velocity_urto_meno=velocities[-1][-1] if psoln is not None else 1
    velocity_urto_piu=r*velocity_urto_meno

    y_curr=[0, velocity_urto_piu]
    if activation>0:
        ts.append([])
        ts[-1] = np.arange(ts[-2][len(rotations[-1])], tStop, tInc)
        interpol=interpolate.interp1d(t[len(rotations[-1]):],acc[len(rotations[-1]):],fill_value="extrapolate")
    else:
        ts.append([])
        ts[-1] = np.arange(0. + ndv*tInc, tStop, tInc)
        interpol=interpolate.interp1d(t[ndv:],acc[ndv:],fill_value="extrapolate")
        #interpol=interpolate.interp1d(t[len(rotations[-1]):],acc[len(rotations[-1]):],fill_value="extrapolate")
    params = [p_quadro, alpha, interpol,g]
    psoln = odeint(f, y_curr, ts[-1], args=(params,),full_output=True) 

#%%
    
numero_di_volte=1
for ndv in range(numero_di_volte):
    rotations.append([])
    velocities.append([])
    
    if psoln is not None:
        for i in range(len(psoln[0][:,0])):
            if psoln[0][i,0]>=0:
                rotations[-1]=rotations[-1]+[psoln[0][i,0]*180/np.pi]
                velocities[-1]=velocities[-1]+[psoln[0][i,1]]
                #ndv_list.append(ndv) if ndv not in ndv_list else 0
            else:
                break   
    
    r=e_1s if velocities[-1] else 0
    velocity_urto_meno=velocities[-1][-1] if velocities[-1] else 1
    velocity_urto_piu=r*velocity_urto_meno

    ts.append([])
    ts[-1]= np.linspace(ts[-2][len(rotations[-1])], tStop, int(tStop/tInc)+1)
    y_curr=[0, velocity_urto_piu]
    
    interpol=interpolate.interp1d(t[len(rotations[-1]):],acc[len(rotations[-1]):],fill_value="extrapolate")
    print(ts[-1][0],acc[len(rotations[-1])],y_curr)
    params = [p_quadro, alpha, interpol,g]
    psoln = odeint(f, y_curr, ts[-1], args=(params,),full_output=True) 


#%% tentativi precedenti

#ground motion
#theta in function of time for accelerogram ground motion
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import xlrd
from scipy import interpolate 
    
b=0.3 #mezza base
h=2.4 #mezza altezza
H=2*h
zeta=h/b
vol=2*b*2*h
density=203
g=9.81
m=vol*density
R=np.sqrt((b)**2+(h)**2)
I=(4/3)*m*R**2
p_quadro=m*g*R/I
p=np.sqrt(p_quadro)
alpha=np.arctan(b/h)

# Make time array for solution
tStop = 8.000
tInc = 0.005
t = np.linspace(0., tStop, int(tStop/tInc)+1)

#seismic action
wb=xlrd.open_workbook('C:\\Users\\victor\\Desktop\\PHD\\rocking\\acceleration_laquila.xlsx')
s=wb.sheet_by_name('acceleration strong')
s.nrows
#time=np.linspace(0., tStop, int(tStop/tInc)+1)
#time=time[6000:10000]
acc=[s.cell_value(i,4)*0.01 for i in range(s.nrows)]

# we assume that time and acc are of the same dimension
#dict_time_acc = {}
#for i in range(len(time)):
#    dict_time_acc[time[i]] = acc[i]
#
#def interpolate_f(t):
#
#    return dict_time_acc[t]
    

t=t[400:]
#acc=acc[6000:10000]
#extrapolate creates values after the last one
interpol=interpolate.interp1d(t,acc,fill_value="extrapolate") 
#plt.plot(time,acc)
#plt.plot(t,interpol(t))
#plt.plot(t,[interpolate_f(x) for x in t])
plt.plot(t,interpol(t))

fig, ax = plt.subplots()
axes = plt.gca()
ax.plot(t,interpol(t))
ax.set(xlabel='Time [m]',ylabel='acc [m/s^2]')
#axes.set_ylim([-10,10]) 
axes.set_xlim([0,6]) 

theta_plot = []
t_plot=[]

def f(y, t, params): 
    theta,omega = y
    #print(theta)
    p_quadro, alpha, interpol,g = params
    if theta<0:
        break
    t_plot.append(t)
    theta_plot.append(theta)
    
    derivs = [omega,
              -p_quadro*(np.sin(alpha-theta))+10*p_quadro*interpol(t)*np.cos(alpha-theta)/g]
    return derivs

# Initial values
theta0 = 0*np.pi/180     # initial angular displacement RADIANTI
omega0 = 0.0  # initial angular velocity
params = [p_quadro, alpha, interpol,g]
y0 = [theta0, omega0]
# Call the ODE solver
psoln = odeint(f, y0, t, args=(params,), full_output=True)

# Plot theta as a function of time
theta_plot=np.array(theta_plot)
t_plot=np.array(t_plot)
fig, ax = plt.subplots()
axes = plt.gca()
#ax.plot(t, psoln[0][:,0]*180/np.pi)
ax.plot(t,interpol(t))
ax.plot(t_plot, theta_plot*180/np.pi)
ax.set(xlabel='Time [m]',ylabel='rotations [degree]')
#axes.set_ylim([-5,2]) 
#axes.set_xlim([0,6]) 
#%%
#LUCA
def f(y, t, params):
    theta, omega = copy.deepcopy(y)
    print(y)
    print(omega)
    print(t)
    p_quadro, alpha = params
#    if theta<0:
#        theta = 0
#        omega = e_1s*omega
#        print(theta)
#        print(omega)
#        print()

#    print(theta)
#    print(omega)
#    print()
    derivs = [omega,
              -p_quadro*(np.sin(alpha-theta))]#+10*p_quadro*interpol(t)*np.cos(alpha-theta)/g]    
    return derivs

params = [p_quadro, alpha]


