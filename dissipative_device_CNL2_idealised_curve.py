import numpy as np
import matplotlib.pyplot as plt
from numpy import trapz

B=0.6
H=4.8
R=np.sqrt((B/2)**2+(H/2)**2)
l_h=1
A=B*l_h
ro=203 #massa Kg/m^3
m,g=(B*H*l_h*ro,9.8)
#m,g,R= (sym.Symbol('m'),sym.Symbol('g'),sym.Symbol('R'))
alfa= np.arctan(B/H) #circa 7 gradi
teta_vector=np.arange(0,alfa,0.0002)

#u_teta
k_n=6e6+6e6*0.5 #Pa
#f_m=2*m*g/(0.35*l_h)
f_m=32739+32739*0.5 #Pa
teta_j0=2*m*g/(k_n*B**2*l_h)
a_c=2*m*g/(f_m*l_h)
teta_jc=f_m/(k_n*a_c)
def u_teta(teta_vector):           
    u_teta_dict={}
    teta_c_val=[]
    teta_c_index=[]
    for teta in teta_vector:
        if teta<teta_j0:
            u_teta=B/2-(B**3*k_n*l_h*teta/(12*m*g))
        if teta <=teta_jc and teta >teta_j0: 
            u_teta=(1/3)*np.sqrt(2*m*g/(k_n*l_h*teta))
        if teta>teta_jc:
            if len(teta_c_val)==0:
                teta_c_val.append(teta)
            if len(teta_c_index)==0:
                teta_c_index.append(list(teta_vector).index(teta))
            u_teta=0.5*(m*g/(f_m*l_h)+(f_m**3*l_h)/(12*m*g*k_n**2*teta**2))
        u_teta_dict[teta]=u_teta
    return u_teta_dict, teta_c_val[0], teta_c_index[0]
            
u_teta_dict=u_teta(teta_vector)[0]
teta_c_val=u_teta(teta_vector)[1]
teta_c_index=u_teta(teta_vector)[2]
u_teta=np.array([k for k in u_teta_dict.values()])        
plt.plot(teta_vector,u_teta)

def gamma_incognita(teta_vector,u_teta):
    gamma_dict_flex={}
    gamma_dict_rigid={}
    for i, teta in enumerate(teta_vector):
        M_stab_flex=0
        M_stab_rigid=0
        M_stab_flex=m*g*R*np.sin(alfa-teta)-m*g*u_teta[i]
        M_stab_rigid=m*g*R*np.sin(alfa-teta)
        M_ovt=m*g*R*np.cos(alfa-teta)
        gamma_flex=M_stab_flex/M_ovt
        gamma_rigid=M_stab_rigid/M_ovt
        gamma_dict_flex[teta]=gamma_flex
        gamma_dict_rigid[teta]=gamma_rigid
    return gamma_dict_flex, gamma_dict_rigid

gamma_dict_flex=gamma_incognita(teta_vector,u_teta)[0]
gamma_dict_rigid=gamma_incognita(teta_vector,u_teta)[1]
gamma_flex=np.array([k for k in gamma_dict_flex.values()])    
gamma_rigid=np.array([k for k in gamma_dict_rigid.values()])    
plt.plot(teta_vector,gamma_flex)
plt.plot(teta_vector,gamma_rigid)
#%%
#anchor
L_t=4
H_t=H-0.3
epsilon_y=0.002 #0.2%
epsilon_ul=0.016 #10%
sigma_y=235e3 #KN/m^2
E=sigma_y/epsilon_y # 117.5e6 KN/m^2
fi_anchor=0.015 #m
num_anchors=4
A=num_anchors*np.pi*fi_anchor**2/4 #m^2
actual_force=sigma_y*A
teta_y=np.arctan(epsilon_y*L_t/H_t)
teta_ul=np.arctan(epsilon_ul*L_t/H_t)

#introduco il device
F_sliding=actual_force*0.90
u_sliding=0.04
epsilon_dev1=F_sliding/(A*E)
epsilon_dev2=epsilon_dev1+u_sliding/L_t
sigma_dev=F_sliding/A
epsilon_ul_dev=u_sliding/L_t+epsilon_ul
epsilon_y_dev=u_sliding/L_t+epsilon_y
#sigma_y
#epsilon_y
#epsilon_ul

tratto1=np.arange(0,epsilon_dev1,0.0001)
tratto2=np.arange(tratto1[-1],epsilon_dev2,0.0001)
tratto3=np.arange(tratto2[-1],u_sliding/L_t+epsilon_y,0.0001)
tratto4=np.arange(tratto3[-1],u_sliding/L_t+epsilon_ul,0.0001)
concateno_tratti_epsilon=np.concatenate((tratto1,tratto2,tratto3,tratto4))

sigma_tratto1=E*tratto1
sigma_tratto2=np.ones(len(tratto2))*E*epsilon_dev1
sigma_tratto3=sigma_tratto2[-1]+E*(tratto3-tratto2[-1])
sigma_tratto4=np.ones(len(tratto4))*E*epsilon_y
concateno_tratti_sigma=np.concatenate((sigma_tratto1,sigma_tratto2,sigma_tratto3,sigma_tratto4))

plt.plot(concateno_tratti_epsilon,concateno_tratti_sigma)

fig, ax = plt.subplots()
axes = plt.gca()
ax.plot(concateno_tratti_epsilon,concateno_tratti_sigma*A, label='Inelastic demand A-D spectrum',linewidth=2, color='black')
plt.xticks([epsilon_dev1, epsilon_dev2,u_sliding/L_t+epsilon_y+0.001,u_sliding/L_t+epsilon_ul],('$\\theta_{d1}$','$\\theta_{d2}$',\
           '$\\theta_y$', '$\\theta_{ul}$'))
plt.yticks([sigma_dev*A, sigma_y*A],('$F_{sliding}$','$F_{yield}$'))
ax.tick_params(axis='both', which='major', labelsize=12)
#ax.set(xlabel='$\epsilon $',ylabel='F')
axes.set_ylim([0,sigma_y*A+5]) 
axes.set_xlim([0,u_sliding/L_t+epsilon_ul+0.003]) 
fig.savefig("output_python_figure\constitutive law device and anchor_teta.png", dpi = 500, bbox_inches='tight')
#%%

#introduction of device in equation
teta_dev1=np.arctan(epsilon_dev1*L_t/H_t)
teta_dev2=np.arctan(epsilon_dev2*L_t/H_t)
teta_y_dev=np.arctan((u_sliding/L_t+epsilon_y)*L_t/H_t)
teta_ul_dev=np.arctan((u_sliding/L_t+epsilon_ul)*L_t/H_t)

def F_dev(teta_vector):
    F_dev_dict={}
    teta_dev1_val=[]
    teta_dev1_index=[]
    teta_dev2_val=[]
    teta_dev2_index=[]
    teta_yield_dev=[]
    teta_yield_dev_index=[]
    teta_plastic_dev=[]
    teta_plastic_dev_index=[]
    for i, teta in enumerate(teta_vector):
        if teta<teta_dev1:
            F_dev=H_t*np.tan(teta)/L_t*E*A
        if teta>teta_dev1 and teta<teta_dev2:
            F_dev=E*epsilon_dev1*A
            if len(teta_dev1_val)==0:
                teta_dev1_val.append(teta)
            if len(teta_dev1_index)==0:
                teta_dev1_index.append(list(teta_vector).index(teta))
        if teta>=teta_dev2 and teta<teta_y_dev:
            F_dev=E*epsilon_dev1*A+H_t*np.tan(teta-teta_dev2)/L_t*E*A
            if len(teta_dev2_val)==0:
                teta_dev2_val.append(teta)
            if len(teta_dev2_index)==0:
                teta_dev2_index.append(list(teta_vector).index(teta))
        if teta>teta_y_dev and teta<=teta_ul_dev:
            F_dev=E*epsilon_y*A
            if len(teta_yield_dev)==0:
                teta_yield_dev.append(teta)
            if len(teta_yield_dev_index)==0:
                teta_yield_dev_index.append(list(teta_vector).index(teta))
        if teta>teta_ul_dev:
            F_dev=0
            if len(teta_plastic_dev)==0:
                teta_plastic_dev.append(teta)
            if len(teta_plastic_dev_index)==0:
                teta_plastic_dev_index.append(list(teta_vector).index(teta)-1)
        F_dev_dict[teta]=F_dev
    return F_dev_dict, teta_dev1_val[0],teta_dev1_index[0],teta_dev2_val[0],teta_dev2_index[0],teta_yield_dev[0],teta_yield_dev_index[0],teta_plastic_dev[0], teta_plastic_dev_index[0]

F_dev_dict,teta_dev1_val,teta_dev1_index,teta_dev2_val,teta_dev2_index,teta_yield_dev,teta_yield_dev_index,teta_plastic_dev, teta_plastic_dev_index=F_dev(teta_vector)
F_dev=np.array([k for k in F_dev_dict.values()]) 
plt.plot(teta_vector,F_dev)
     
def gamma_incognita_device(teta_vector,u_teta,F_dev):
    gamma_dict_flex_dev={}
    gamma_dict_rigid_dev={}
    for i, teta in enumerate(teta_vector):
        M_stab_flex_dev=0
        M_stab_rigid_dev=0
        M_stab_flex_dev=m*g*R*np.sin(alfa-teta)-m*g*u_teta[i]+F_dev[i]*H_t
        M_stab_rigid_dev=m*g*R*np.sin(alfa-teta)+F_dev[i]*H_t
        M_ovt=m*g*R*np.cos(alfa-teta)
        gamma_flex_dev=M_stab_flex_dev/M_ovt
        gamma_rigid_dev=M_stab_rigid_dev/M_ovt
        gamma_dict_flex_dev[teta]=gamma_flex_dev
        gamma_dict_rigid_dev[teta]=gamma_rigid_dev
    return gamma_dict_flex_dev, gamma_dict_rigid_dev

gamma_dict_flex_dev=gamma_incognita_device(teta_vector,u_teta,F_dev)[0]
gamma_dict_rigid_dev=gamma_incognita_device(teta_vector,u_teta,F_dev)[1]
gamma_flex_dev=np.array([k for k in gamma_dict_flex_dev.values()])    
gamma_rigid_dev=np.array([k for k in gamma_dict_rigid_dev.values()])    

#punto di controllo a mezzeria
fig, ax = plt.subplots()
axes = plt.gca()
ax.plot(H/2*np.tan(teta_vector),gamma_flex_dev/gamma_rigid[0], label='flexible interface + device',linewidth=1, color='blue')
#ax.plot(H/2*np.tan(teta_vector),gamma_flex_dev/gamma_rigid_dev[0], label='flexible interface',linewidth=1, color='red')
ax.plot(H/2*np.tan(teta_vector),gamma_rigid/gamma_rigid[0], label='rigid interface',linewidth=1, color='black',linestyle='--')
ax.set(xlabel='Rotation $\\theta $ ',ylabel='$\lambda$')
axes.set_ylim([0,1]) 
axes.set_xlim([0,0.4]) 
legend = ax.legend(loc='upper right', shadow=False, fontsize='medium')
#fig.savefig("output_python_figure\capacity curves.png", dpi = 500, bbox_inches='tight')

#in accelerazioni e displacements
#spostamento del punto di controllo (device)

#NON NORMALIZZATO
displacements_CP=H_t*np.tan(teta_vector)
fig, ax = plt.subplots()
axes = plt.gca()
ax.plot(displacements_CP,gamma_flex_dev/gamma_rigid[0], label='flexible interface + tie-rod',linewidth=1, color='blue')
#ax.plot(displacements_CP,gamma_flex/gamma_rigid[0], label='flexible interface',linewidth=1, color='red')
ax.plot(displacements_CP,gamma_rigid/gamma_rigid[0], label='rigid interface',linewidth=1, color='black',linestyle='--')
ax.plot([H_t*np.tan(teta_c_val), H_t*np.tan(teta_c_val)], [0,gamma_flex_dev[teta_c_index]/gamma_rigid[0]], label='rigid interface',linewidth=1, color='black',linestyle='--')
ax.set(xlabel='Displacements [m] ',ylabel='$\lambda$')
axes.set_ylim([0,1]) 
axes.set_xlim([0,0.7]) 
legend = ax.legend(loc='upper right', shadow=False, fontsize='medium')
#fig.savefig("output_python_figure\capacity curve in diplacemetn of anchor.png", dpi = 500, bbox_inches='tight')
#%%
#spettro per il design
#elastic response spectrum
a_g,eta,g,F_0=(0.24,1,9.81,2.2)

S=1.15

#soil Type B
Tb=0.2
Tc=0.6
Td=2.00
seconds=4
T=np.arange(seconds,step=seconds/399) #len(d_star_anchoring)


def accelerazione_spettrale(T):
    S_e=[]
    for t in T:
        if t <=Tb:
            S_e+=[a_g*S*(1+t/Tb*(eta*F_0-1))]
        if t <=Tc and t >Tb: # tratto orizontale
            S_e+=[a_g*S*eta*F_0]
        if t >Tc and t <Td:
            S_e+=[a_g*S*eta*F_0*Tc/t]
        if t >Td:
            S_e+=[a_g*S*eta*F_0*Tc*Td/t**2]
    S_e=np.array(S_e)
    return S_e
S_e=accelerazione_spettrale(T)    

#valori frnacesco morfuni
#T_el_FM=1.15 #s
#S_e_FM=a_g*S*(1+T_el_FM/Tb*(eta*F_0-1))

fig, ax = plt.subplots()
axes = plt.gca()
ax.plot(T,S_e*g, label='Elastic spectrum',linewidth=1, color='black')
ax.set(xlabel='Period [s]',ylabel='Spectral acceleration [m/s^2]')
#axes.set_ylim([0,0.8]) 
#axes.set_xlim([0,3.7]) 
#ax.xaxis.set_ticks(np.arange(0,3.7,0.3))
#ax.yaxis.set_ticks(np.arange(0,8.1,2))
legend = ax.legend(loc='upper right', shadow=False, fontsize='medium')
#fig.savefig("output_python_figure\design spectra.png", dpi = 500, bbox_inches='tight')

def spostamento_spettrale(S_e):
    S_de=[S_e[i]*g*(T[i]/(2*np.pi))**2 for i in range(len(S_e))]
    return np.array(S_de)

S_de=spostamento_spettrale(S_e)

fig, ax = plt.subplots()
axes = plt.gca()
ax.plot(S_de,S_e*g, label='Elastic demand A-D spectrum',linewidth=1, color='black')
ax.set(xlabel='Spectral displacement [m]',ylabel='Spectral acceleration [m/s^2]')
axes.set_ylim([0,9]) 
axes.set_xlim([-0.002,0.3]) 
ax.xaxis.set_ticks(np.arange(0,0.3,0.05))
ax.yaxis.set_ticks(np.arange(0,9.1,2))
legend = ax.legend(loc='upper right', shadow=False, fontsize='medium')
#fig.savefig("output_python_figure\Elastic demand A-D spectrum.png", dpi = 500, bbox_inches='tight')

#%%
S_d_yield_dev=epsilon_dev1*L_t
S_a_yield_dev=gamma_flex_dev[teta_dev1_index]*g
S_d_yield_bar=epsilon_y_dev*L_t
S_a_yield_bar=gamma_flex_dev[teta_yield_dev_index]*g
#%%

#FOCUS OF CAPACITY CURVE AND CRITICAL POINTS
disp_crush_toe=H_t*np.tan(teta_c_val)
gamma_crush_toe=gamma_flex_dev[teta_c_index]
disp_sliding=H_t*np.tan(teta_dev1)
gamma_sliding=gamma_flex_dev[teta_dev1_index]
disp_plastic_fail=H_t*np.tan(teta_ul_dev)
gamma_ultimate=gamma_flex_dev[teta_plastic_dev_index]

fig, ax = plt.subplots()
axes = plt.gca()
ax.plot(displacements_CP[:teta_plastic_dev_index+1],gamma_flex_dev[:teta_plastic_dev_index+1]*g, label="System's capacity",linewidth=1, color='blue')
ax.plot(displacements_CP[teta_plastic_dev_index:],gamma_flex[teta_plastic_dev_index:]*g, linewidth=0.6, color='blue',linestyle='--')

ax.plot([disp_sliding, disp_sliding], [0,gamma_sliding*g],linewidth=1, color='black',linestyle='--')
ax.plot([S_d_yield_bar,S_d_yield_bar],[0,S_a_yield_bar],linewidth=1, color='black',linestyle='--')
#ax.plot([disp_crush_toe, disp_crush_toe], [0,gamma_crush_toe*g],linewidth=1, color='black',linestyle='--')
ax.plot([disp_plastic_fail, disp_plastic_fail], [0,gamma_ultimate*g],linewidth=1, color='black',linestyle='--')
size=25
ax.scatter(disp_sliding,gamma_sliding*g,linewidth=1, color='blue', s=size)
ax.scatter(S_d_yield_bar,S_a_yield_bar,linewidth=1, color='blue', s=size)
ax.scatter(disp_crush_toe,gamma_crush_toe*g,linewidth=1, color='blue', s=size)
ax.scatter(disp_plastic_fail,gamma_ultimate*g,linewidth=1, color='blue', s=size)
plt.annotate('Toe\n crush',xy=(disp_crush_toe,gamma_crush_toe*g-0.18), xytext=(disp_crush_toe,gamma_crush_toe*g-0.18),fontsize=10)

plt.xticks([disp_sliding, S_d_yield_bar,disp_plastic_fail],('Sliding','Yielding', 'Tie\n break'),horizontalalignment='center')
ax.set(xlabel='$ S_d (T)$ [m]',ylabel='$ S_a$(T) [$m/s^2$]')
axes.set_ylim([0,gamma_ultimate*g+0.3]) 
axes.set_xlim([0,0.12]) 
#ax.xaxis.set_ticks(np.arange(0,0.161,0.1))
#ax.yaxis.set_ticks(np.arange(0,9,4))
legend = ax.legend(loc='upper right', shadow=False, fontsize='medium') 
#%%
#NUMERICAL INTERGATION AND LINEARIZATION TO ELASTO PERFECTLY PLASTIC


if disp_plastic_fail<disp_crush_toe:
    disp_critical_ropture=disp_plastic_fail
    acc_critical_ropture=gamma_ultimate
    critical_ropture_index=teta_plastic_dev_index-1
else:
    disp_critical_ropture=disp_crush_toe
    acc_critical_ropture=gamma_crush_toe
    critical_ropture_index=teta_c_index
    
x=displacements_CP[:critical_ropture_index]
y=gamma_flex_dev[:critical_ropture_index]*g
# Compute the area using the composite trapezoidal rule.
area = trapz(y,x) #integra sull asse x
print("area =", area)
#elasto plastic values
delta_y=2*(disp_critical_ropture-area/(acc_critical_ropture*g))
S_a_y=acc_critical_ropture*g

fig, ax = plt.subplots()
axes = plt.gca()
ax.plot(displacements_CP[:critical_ropture_index+1],gamma_flex_dev[:critical_ropture_index+1]*g, label="System's capacity ",linewidth=1, color='blue')
ax.plot(displacements_CP[critical_ropture_index:],gamma_flex_dev[critical_ropture_index:]*g, linewidth=0.8, color='blue',linestyle='--',label="System's capacity after failure ")
ax.plot([0,delta_y,disp_critical_ropture],[0,S_a_y,S_a_y],linewidth=2, color='red',label="Idealised Curve")
        
ax.plot([disp_sliding, disp_sliding], [0,gamma_sliding*g],linewidth=1, color='black',linestyle='--')
ax.plot([S_d_yield_bar,S_d_yield_bar],[0,S_a_yield_bar],linewidth=1, color='black',linestyle='--')
ax.plot([disp_crush_toe, disp_crush_toe], [0,gamma_crush_toe*g],linewidth=1, color='black',linestyle='--')
ax.plot([disp_plastic_fail, disp_plastic_fail], [0,gamma_ultimate*g],linewidth=1, color='black',linestyle='--')
size=25
ax.scatter(disp_sliding,gamma_sliding*g,linewidth=1, color='blue', s=size)
ax.scatter(S_d_yield_bar,S_a_yield_bar,linewidth=1, color='blue', s=size)
ax.scatter(disp_crush_toe,gamma_crush_toe*g,linewidth=1, color='blue', s=size+2)
ax.scatter(disp_plastic_fail,gamma_ultimate*g,linewidth=1, color='red', s=size)

ax.scatter(delta_y,S_a_y,linewidth=1, color='red', s=size+2)
plt.annotate('$\Delta_y$', xy=(delta_y,S_a_y+0.15), xytext=(delta_y,S_a_y+0.025),fontsize=13)
plt.annotate('$\Delta_u$', xy=(disp_crush_toe,S_a_y+0.15), xytext=(disp_crush_toe,S_a_y+0.025),fontsize=13)

plt.annotate('Toe\n crush',xy=(disp_crush_toe,gamma_crush_toe*g-0.18), xytext=(disp_crush_toe,gamma_crush_toe*g-0.2),fontsize=10)

plt.xticks([disp_sliding, S_d_yield_bar,disp_plastic_fail],('Sliding','Yielding', 'Tie\n break'),horizontalalignment='center')
ax.set(xlabel='$ S_d (T)$ [m]',ylabel='$ S_a$(T) [$m/s^2$]')
axes.set_ylim([0,gamma_ultimate*g+0.6]) 
axes.set_xlim([0,0.12]) 
#ax.xaxis.set_ticks(np.arange(0,0.161,0.1))
#ax.yaxis.set_ticks(np.arange(0,9,4))
#ax.xaxis.set_ticks(np.arange(0,0.161,0.1))
#ax.yaxis.set_ticks(np.arange(0,9,4))
legend = ax.legend(loc='upper right', shadow=False, fontsize='medium') 
#fig.savefig("output_python_figure\idealised curve.png", dpi = 500, bbox_inches='tight')
#%%
#plot together the demand and the capacity
#transformation into spectral values
mu_capacity=disp_critical_ropture/delta_y
S_a_yield=acc_critical_ropture*g
T_el=2*np.pi*np.sqrt(delta_y/S_a_yield)

#spostamenti virtuali
W=m*g
teta_virt=1/H_t
delta_x_T_anchor=1 
delta_x_W= H/2*teta_virt
delta_x=[delta_x_W]
loads=[W]
#Massa partecipante al cinematismo
def M_star(loads, delta_x):
    numeratore=0
    denominatore=0
    for i in range(len(loads)):
        numeratore+=loads[i]*delta_x[i]
        denominatore+= loads[i]*((delta_x[i])**2)
    return numeratore**2/(g*denominatore)
M_star=M_star(loads, delta_x) #KN

#Frazione di massa partecipante
def e_star(loads, M_star):
    denominatore=0
    for i in range(len(loads)):
        denominatore+= loads[i]
    return g*M_star/denominatore
e_star=e_star(loads, M_star)
S_d_PP_dev=disp_critical_ropture

S_d_yield_dev=delta_y
S_a_yield_dev=acc_critical_ropture*g

def S_a_el(T_el):
    S_a_el=0
    t=T_el
    if t <=Tb:
        S_a_el=a_g*S*(1+t/Tb*(eta*F_0-1))
    if t <=Tc and t >Tb: # tratto orizontale
        S_a_el=a_g*S*eta*F_0
    if t >Tc and t <Td:
        S_a_el=a_g*S*eta*F_0*Tc/t
    if t >Td:
        S_a_el=a_g*S*eta*F_0*Tc*Td/t**2
    return S_a_el*g
S_a_el_y=S_a_el(T_el) #4.24
S_d_el_y=S_a_el_y*(T_el/(2*np.pi))**2 #0.15

fig, ax = plt.subplots()
axes = plt.gca()
ax.plot(S_de,S_e*g, label='Elastic demand A-D spectrum',linewidth=1, color='black')
ax.plot(displacements_CP[:teta_plastic_dev_index+2],gamma_flex_dev[:teta_plastic_dev_index+2]*g, label='flexible interface + tie-rod',linewidth=1, color='blue')
ax.plot(displacements_CP[teta_plastic_dev_index:],gamma_flex[teta_plastic_dev_index:]*g, linewidth=1, color='blue')
ax.plot([H_t*np.tan(teta_c_val), H_t*np.tan(teta_c_val)], [0,gamma_flex_dev[teta_c_index]/gamma_rigid[0]], label='rigid interface',linewidth=1, color='black',linestyle='--')
ax.plot([0,delta_y,disp_critical_ropture],[0,S_a_y,S_a_y],linewidth=2, color='red',label="Idealised Curve")

ax.plot([0,S_d_el_y],[0,S_a_el_y], color='red',linestyle='--')

ax.set(xlabel='$ S_d (T)$ [m]',ylabel='$ S_a$(T) [$m/s^2$]')
axes.set_ylim([0,9]) 
axes.set_xlim([-0.002,0.41]) 
ax.xaxis.set_ticks(np.arange(0,0.41,0.1))
ax.yaxis.set_ticks(np.arange(0,9.1,2))
legend = ax.legend(loc='upper right', shadow=False, fontsize='medium')      
#fig.savefig("output_python_figure\Elastic demand A-D spectrum_with_capacity curves.png", dpi = 500, bbox_inches='tight')
#%% 

mu_demand=mu_capacity
def accelerazione_spettrale_inelastica(T):
    S_e_inel=[]

    for t in T:
        if t<Tc:
            R_mu= ((mu_demand-1)*t/Tc)+1
        if t>=Tc:
            R_mu= mu_demand
        if t <=Tb:
            S_e_inel+=[a_g*S*(1+t/Tb*(eta*F_0-1))/R_mu]
        if t <=Tc and t >Tb: # tratto orizontale
            S_e_inel+=[a_g*S*eta*F_0/R_mu]
        if t >Tc and t <Td:
            S_e_inel+=[a_g*S*eta*F_0*Tc/t/R_mu]
        if t >Td:
            S_e_inel+=[a_g*S*eta*F_0*Tc*Td/t**2/R_mu]
    S_e_inel=np.array(S_e_inel)
    return S_e_inel, R_mu

S_e_inel,R_mu=accelerazione_spettrale_inelastica(T)     

def spostamento_spettrale_inel(S_de,T):
    S_de_inel=[]  
    for i, t in enumerate(T):
        if t<Tc:
            R_mu= ((mu_demand-1)*t/Tc)+1
        if t>=Tc:
            R_mu= mu_demand
        S_de_inel+=[S_de[i]*mu_demand/R_mu]
    return np.array(S_de_inel)

S_de_inel=spostamento_spettrale_inel(S_de,T)   

fig, ax = plt.subplots()
axes = plt.gca()
ax.plot(S_de_inel,S_e_inel*g, label='Inelastic demand A-D spectrum',linewidth=2, color='black')
ax.plot(S_de,S_e*g, label='Elastic demand A-D spectrum',linewidth=1, color='brown')
ax.plot(displacements_CP[:teta_plastic_dev_index+2],gamma_flex_dev[:teta_plastic_dev_index+2]*g, label='flexible interface + tie-rod',linewidth=1, color='blue')
ax.plot(displacements_CP[teta_plastic_dev_index:],gamma_flex[teta_plastic_dev_index:]*g, linewidth=0.6, color='blue',linestyle='--')
ax.plot([H_t*np.tan(teta_c_val), H_t*np.tan(teta_c_val)], [0,gamma_flex_dev[teta_c_index]*g],linewidth=1, color='black',linestyle='--')
ax.plot([0,delta_y,disp_critical_ropture],[0,S_a_y,S_a_y],linewidth=2, color='red',label="Idealised Curve")

ax.plot([0,S_d_el_y],[0,S_a_el_y], color='red',linestyle='--') #elastic system
#ax.plot([S_d_el_y,S_d_el_y],[0,S_a_el_y],linewidth=1, color='black',linestyle='--')
ax.plot([S_d_yield_bar,S_d_yield_bar],[0,S_a_yield_bar],linewidth=1, color='black',linestyle='--')
ax.set(xlabel='$ S_d (T)$ [m]',ylabel='$ S_a$(T) [$m/s^2$]')
axes.set_ylim([0,7]) 
axes.set_xlim([-0.002,0.21]) 
ax.xaxis.set_ticks(np.arange(0,0.21,0.1))
ax.yaxis.set_ticks(np.arange(0,7.1,2))
legend = ax.legend(loc='upper right', shadow=False, fontsize='medium')     
plt.annotate('$R_{\mu}$ = '+'{:.1f}'.format(mu_capacity),xy=(0.05,3), xytext=(0.022,3),fontsize=11)
plt.annotate('Toe crush',xy=(disp_crush_toe+0.003,0.1), xytext=(disp_crush_toe+0.003,0.1),fontsize=11)
plt.annotate('Yielding',xy=(S_d_yield_bar+0.003,0.1), xytext=(S_d_yield_bar+0.003,0.1),fontsize=11)

#fig.savefig("output_python_figure\inelastic demand A-D spectrum_with_capacity curves_device.png", dpi = 500, bbox_inches='tight')



 