import numpy as np
import matplotlib.pyplot as plt
from numpy import trapz
#%%
#FAIL PARAMETERS
k_n=6e6+6e6*0.5 #Pa
f_m=32739+32739*0.5 #Pa 
fi_anchor=0.015+0.015*0 #m
#%%
#RETROFIT PARAMETERS
k_n=6e6+6e6*0.7 #Pa
f_m=32739+32739*0.7 #Pa
fi_anchor=0.015+0.015*0.3 #m
#%%

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

#k_n=6e6+6e6*0.5 #Pa
#f_m=32739+32739*0.5 #Pa
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
#introduco l'anchor
L_t=4
H_t=H-0.3
epsilon_y=0.002 #0.2%
epsilon_ul=0.016 #10%
sigma_y=235e3 #KN/m^2
E=sigma_y/epsilon_y # 117.5e6 KN/m^2
#fi_anchor=0.015+0.015*0 #m
num_anchors=4
A=num_anchors*np.pi*fi_anchor**2/4 #m^2
sigma_y*A
teta_y=np.arctan(epsilon_y*L_t/H_t)
teta_ul=np.arctan(epsilon_ul*L_t/H_t)

def F_t(teta_vector):
    F_t_dict={}
    teta_yield=[]
    teta_yield_index=[]
    teta_plastic=[]
    teta_plastic_index=[]
    for i, teta in enumerate(teta_vector):
        if teta<teta_y:
            F_t=(H_t*np.tan(teta)/L_t)*E*A
        if teta>=teta_y and teta<teta_ul:
            F_t=epsilon_y*E*A #KN
            if len(teta_yield)==0:
                teta_yield.append(teta)
            if len(teta_yield_index)==0:
                teta_yield_index.append(list(teta_vector).index(teta))
        if teta>teta_ul:
            F_t=0
            if len(teta_plastic)==0:
                teta_plastic.append(teta)
            if len(teta_plastic_index)==0:
                teta_plastic_index.append(list(teta_vector).index(teta))
        F_t_dict[teta]=F_t
    return F_t_dict, teta_yield[0], teta_yield_index[0],teta_plastic[0], teta_plastic_index[0]

F_t_dict,teta_yield,teta_yield_index,teta_plastic,teta_plastic_index=F_t(teta_vector)
F_t=np.array([k for k in F_t_dict.values()]) 
plt.plot(teta_vector,F_t)
     
def gamma_incognita_anchor(teta_vector,u_teta,F_t):
    gamma_dict_flex={}
    gamma_dict_rigid={}
    for i, teta in enumerate(teta_vector):
        M_stab_flex=0
        M_stab_rigid=0
        M_stab_flex=m*g*R*np.sin(alfa-teta)-m*g*u_teta[i]+F_t[i]*H_t
        M_stab_rigid=m*g*R*np.sin(alfa-teta)+F_t[i]*H_t
        M_ovt=m*g*R*np.cos(alfa-teta)
        gamma_flex=M_stab_flex/M_ovt
        gamma_rigid=M_stab_rigid/M_ovt
        gamma_dict_flex[teta]=gamma_flex
        gamma_dict_rigid[teta]=gamma_rigid
    return gamma_dict_flex, gamma_dict_rigid

gamma_dict_flex_anchor=gamma_incognita_anchor(teta_vector,u_teta,F_t)[0]
gamma_dict_rigid_anchor=gamma_incognita_anchor(teta_vector,u_teta,F_t)[1]
gamma_flex_anchor=np.array([k for k in gamma_dict_flex_anchor.values()])    
gamma_rigid_anchor=np.array([k for k in gamma_dict_rigid_anchor.values()])    

disp_crush_toe=H_t*np.tan(teta_c_val)
acc_crush_toe=gamma_flex_anchor[teta_c_index]

fig, ax = plt.subplots()
axes = plt.gca()
ax.plot(H/2*np.tan(teta_vector),gamma_flex_anchor/gamma_rigid[0], label='flexible interface + tie-rod',linewidth=1, color='blue')
ax.plot(H/2*np.tan(teta_vector),gamma_flex/gamma_rigid[0], label='flexible interface',linewidth=1, color='red')
ax.plot(H/2*np.tan(teta_vector),gamma_rigid/gamma_rigid[0], label='rigid interface',linewidth=1, color='black',linestyle='--')
ax.set(xlabel='Rotation $\\theta $ ',ylabel='$\lambda$')
axes.set_ylim([0,1]) 
axes.set_xlim([0,1]) 
legend = ax.legend(loc='upper right', shadow=False, fontsize='medium')
#fig.savefig("output_python_figure\capacity curves.png", dpi = 500, bbox_inches='tight')

#in accelerazioni e displacements
#spostamento del punto di controllo (anchora)

#NON NORMALIZZATO
displacements_CP=H_t*np.tan(teta_vector)
fig, ax = plt.subplots()
axes = plt.gca()
ax.plot(displacements_CP,gamma_flex_anchor/gamma_rigid[0], label='flexible interface + tie-rod',linewidth=1, color='blue')
ax.plot(displacements_CP,gamma_flex/gamma_rigid[0], label='flexible interface',linewidth=1, color='red')
ax.plot(displacements_CP,gamma_rigid/gamma_rigid[0], label='rigid interface',linewidth=1, color='black',linestyle='--')
ax.set(xlabel='Displacements [m] ',ylabel='$\lambda$')
axes.set_ylim([0,1]) 
axes.set_xlim([0,0.7]) 
legend = ax.legend(loc='upper right', shadow=False, fontsize='medium')
#fig.savefig("output_python_figure\capacity curve in diplacemetn of anchor.png", dpi = 500, bbox_inches='tight')

#%%
#NORMALIZZATO
limite_index = [i for i, x in enumerate(gamma_flex_anchor) if x < 0][1]

fig, ax = plt.subplots()
axes = plt.gca()
ax.plot(displacements_CP[:teta_plastic_index+1]/displacements_CP[-1],gamma_flex_anchor[:teta_plastic_index+1]/gamma_rigid[0], label='flexible interface + tie-rod',linewidth=1, color='blue')
ax.plot(displacements_CP/displacements_CP[-1],gamma_flex/gamma_rigid[0], label='flexible interface',linewidth=1, color='red')
ax.plot(displacements_CP/displacements_CP[-1],gamma_rigid/gamma_rigid[0], label='rigid interface',linewidth=1, color='black',linestyle='--')
ax.plot([disp_crush_toe/displacements_CP[-1], disp_crush_toe/displacements_CP[-1]], [0,acc_crush_toe/gamma_rigid[0]],linewidth=1, color='black',linestyle='--')

ax.set(xlabel='Normalized Drift',ylabel='$\lambda$')
axes.set_ylim([0,1]) 
axes.set_xlim([0,1]) 
legend = ax.legend(loc='upper right', shadow=False, fontsize='medium')
plt.annotate('Toe crush',xy=(disp_crush_toe+0.003,0.2), xytext=(disp_crush_toe/displacements_CP[-1]+0.01,0.2),fontsize=11)

#fig.savefig("output_python_figure\capacity curve of restrained anchor_normalized.png", dpi = 500, bbox_inches='tight')
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
axes.set_ylim([0,7]) 
axes.set_xlim([-0.002,0.21]) 
ax.xaxis.set_ticks(np.arange(0,0.21,0.05))
ax.yaxis.set_ticks(np.arange(0,7,2))
legend = ax.legend(loc='upper right', shadow=False, fontsize='medium')
#fig.savefig("output_python_figure\Elastic demand A-D spectrum.png", dpi = 500, bbox_inches='tight')
#%%
#ZOOM
disp_yielding=displacements_CP[teta_yield_index]
acc_yielding=gamma_flex_anchor[teta_yield_index]
disp_ultimate=displacements_CP[teta_plastic_index-1]
acc_ultimate=gamma_flex_anchor[teta_plastic_index-1]
disp_crush_toe=H_t*np.tan(teta_c_val)
acc_crush_toe=gamma_flex_anchor[teta_c_index]

#compute the area
x=displacements_CP[:teta_plastic_index]
y=gamma_flex_anchor[:teta_plastic_index]*g
# Compute the area using the composite trapezoidal rule.
area = trapz(y,x) #integra sull asse x
print("area =", area)
#elasto plastic values
delta_y=2*(disp_ultimate-area/(acc_ultimate*g))
S_a_y=acc_ultimate*g

fig, ax = plt.subplots()
axes = plt.gca()
ax.plot(displacements_CP[:teta_plastic_index],gamma_flex_anchor[:teta_plastic_index]*g, label="System's capacity",linewidth=1, color='blue')
ax.plot(displacements_CP[teta_plastic_index-1:],gamma_flex[teta_plastic_index-1:]*g, label="System's capacity after failure",linewidth=1, color='blue',linestyle='--')
ax.plot([disp_yielding, disp_yielding], [0,acc_yielding*g],linewidth=1, color='black',linestyle='--')
ax.plot([disp_ultimate,disp_ultimate],[0,acc_ultimate*g],linewidth=1, color='black',linestyle='--')
ax.plot([disp_crush_toe,disp_crush_toe], [0,acc_crush_toe*g],linewidth=1, color='black',linestyle='--')

ax.plot([0,delta_y,disp_ultimate],[0,S_a_y,S_a_y],linewidth=2, color='red',label="Idealised Curve")
size=25
ax.scatter(disp_yielding,acc_yielding*g,linewidth=1, color='blue', s=size)
ax.scatter(disp_ultimate,acc_ultimate*g,linewidth=1, color='red', s=size)
ax.scatter(disp_crush_toe,acc_crush_toe*g,linewidth=1, color='blue', s=size)
ax.scatter(delta_y,S_a_y,linewidth=1, color='red', s=size+2)
plt.annotate('$\Delta_y$', xy=(delta_y,S_a_y+0.15), xytext=(delta_y,S_a_y+0.025),fontsize=13)
plt.annotate('$\Delta_u$', xy=(disp_ultimate,S_a_y+0.15), xytext=(disp_ultimate,S_a_y+0.025),fontsize=13)
plt.xticks([disp_yielding,disp_ultimate, disp_crush_toe],('Yielding', 'Tie\n break','Toe\n crush'),horizontalalignment='center')

ax.set(xlabel='$ S_d (T)$ [m]',ylabel='$ S_a$(T) [$m/s^2$]')
axes.set_ylim([0,(acc_ultimate*g)+0.7]) 
axes.set_xlim([0,disp_crush_toe+0.01]) 
legend = ax.legend(loc='upper right', shadow=False, fontsize='medium')      
#fig.savefig("output_python_figure\idealised curve of tie system.png", dpi = 500, bbox_inches='tight')
#%%
#plot together the demand and the capacity
#transformation into spectral values

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

S_d_PP=epsilon_ul*L_t
S_d_yield=delta_y
S_a_yield=acc_ultimate*g
mu_capacity=S_d_PP/S_d_yield
T_el=2*np.pi*np.sqrt(S_d_yield/S_a_yield)

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
S_a_el=S_a_el(T_el) #4.67
S_d_el=S_a_el*(T_el/(2*np.pi))**2 #0.14

fig, ax = plt.subplots()
axes = plt.gca()
ax.plot(S_de,S_e*g, label='Elastic demand A-D spectrum',linewidth=1, color='black')
ax.plot(displacements_CP[:teta_plastic_index+1],gamma_flex_anchor[:teta_plastic_index+1]*g, label="System's capacity",linewidth=1, color='blue')
ax.plot(displacements_CP[teta_plastic_index:],gamma_flex[teta_plastic_index:]*g,label="System's capacity after failure", linewidth=1, color='blue',linestyle='--')

#ax.plot([H_t*np.tan(teta_c_val), H_t*np.tan(teta_c_val)], [0,gamma_flex_anchor[teta_c_index]*g], label='rigid interface',linewidth=1, color='black',linestyle='--')
ax.plot([0,delta_y,disp_ultimate],[0,S_a_y,S_a_y],linewidth=2, color='red',label="Idealised Curve")

#ax.plot(displacements_CP,gamma_rigid*g, label='rigid interface',linewidth=1, color='black',linestyle='--')
ax.plot([0,S_d_el],[0,S_a_el], color='red',linestyle='--')
ax.set(xlabel='$ S_d (T)$ [m]',ylabel='$ S_a$(T) [$m/s^2$]')
axes.set_ylim([0,7]) 
axes.set_xlim([-0.002,0.21]) 
ax.xaxis.set_ticks(np.arange(0,0.21,0.1))
ax.yaxis.set_ticks(np.arange(0,7.1,2))
legend = ax.legend(loc='upper right', shadow=False, fontsize='medium')    
plt.annotate('$T_{el}$', xy=(S_d_el,S_a_el+0.01), xytext=(S_d_el,S_a_el+0.01),fontsize=13)
 
#fig.savefig("output_python_figure\Elastic demand A-D spectrum_with_capacity curves.png", dpi = 500, bbox_inches='tight')

#%% 
S_d_PP=epsilon_ul*L_t
S_d_yield=delta_y
S_a_yield=acc_ultimate*g
mu_capacity=S_d_PP/S_d_yield

R_mu_init=S_a_el/S_a_yield  #equivalente a S_d_el/S_d_yield

if T_el<Tc:
    mu_demand=((R_mu_init-1)*Tc/T_el)+1
if T_el>Tc:
    mu_demand=R_mu_init
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
ax.plot(displacements_CP[:teta_plastic_index+1],gamma_flex_anchor[:teta_plastic_index+1]*g, label='Capacity curve up to failure',linewidth=1, color='blue')
ax.plot(displacements_CP[teta_plastic_index:],gamma_flex[teta_plastic_index:]*g, linewidth=0.5, color='blue',linestyle='--')

ax.plot([0,delta_y,disp_ultimate],[0,S_a_y,S_a_y],linewidth=2, color='red',label="Idealised Curve")
ax.plot([0,S_d_el],[0,S_a_el], color='red',linestyle='--')
ax.plot([disp_crush_toe, disp_crush_toe], [0,acc_crush_toe*g],linewidth=1, color='black',linestyle='--')

ax.set(xlabel='$ S_d (T)$ [m]',ylabel='$ S_a$(T) [$m/s^2$]')
axes.set_ylim([0,7]) 
axes.set_xlim([-0.002,0.21]) 
ax.xaxis.set_ticks(np.arange(0,0.21,0.1))
ax.yaxis.set_ticks(np.arange(0,7.1,2))
legend = ax.legend(loc='upper right', shadow=False, fontsize='medium')    
plt.annotate('$T_{el}$', xy=(S_d_el,S_a_el+0.01), xytext=(S_d_el,S_a_el+0.01),fontsize=13)
plt.annotate('$R_{\mu}$ = '+'{:.1f}'.format(mu_capacity),xy=(0.06,3), xytext=(0.025,3),fontsize=11)
plt.annotate('Toe crush',xy=(disp_crush_toe+0.003,0.1), xytext=(disp_crush_toe+0.003,0.1),fontsize=11)

#fig.savefig("output_python_figure\inelastic demand A-D spectrum_with_capacity curves_retrofit.png", dpi = 500, bbox_inches='tight')






