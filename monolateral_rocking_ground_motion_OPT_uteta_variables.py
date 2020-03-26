# ground motion
# theta in function of time for accelerogram ground motion
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import xlrd
from scipy import interpolate

def read_accelerations(filename,time_array):
    """
    Function that reads accelerations from a source file

    :param filename: the name of source file
    :return: a list of accelerations
    """
    wb = xlrd.open_workbook(filename)
    s = wb.sheet_by_name("acceleration strong")
    time_array_length= len(time_array)
    accelerations = [s.cell_value(i, 4) * 0.01 for i in range(time_array_length)]

    return accelerations


def interpolate_time_accelerations(time_array, accelerations):
    """
    Function that, given equal length arrays of times and accelerations,
    interpolate the two and generate a continuous function.

    :param time_array:
    :param accelerations:
    :return: interpolation time-accelerations
    """
    assert len(time_array) == len(accelerations)

    # extrapolate creates values after the last one
    time_acceleration_interpolation = interpolate.interp1d(time_array, accelerations, fill_value="extrapolate")

    return time_acceleration_interpolation

def calculate_u_theta(theta):
    """
    function that computes u_theta and p_quadro_theta for each theta obtained by calculate_derivatives
    """
    if theta <= teta_j0:
        u_theta=(B/2-(B**3*k_n*l_h*theta/(12*m*g)))
        
    if teta_j0 < theta <= teta_jc: 
        u_theta=(1/3)*np.sqrt(2*m*g/(k_n*l_h*(theta)))
    
    if theta > teta_jc:
        u_theta=0.5*(m*g/(f_m*l_h)+(f_m**3*l_h)/(12*m*g*k_n**2*(theta)**2))
    
    R_theta=np.sqrt((R*np.cos(alpha))**2+(R*np.sin(alpha)-u_theta)**2)
    I_theta=m/3*R**2+m*R_theta**2
    p_theta_quadro=m*g*R/I_theta
    
    return u_theta, p_theta_quadro
    
def calculate_derivatives(y, t, p_quadro, alpha, time_acceleration_interpolation, gravity):
    """
    Function that scipy.odeint takes in input in order to calculate the
    rotation theta and the angular velocity omega.

    :param y: output of the differential equation
    :param t: input of the differential equation: array of equal-spaced float
    :param p_quadro: = mgR/I (m = mass of the wall, g: gravity acceleration, R = distance between the centroid G
                              and the rotation point O, I: polar moment of inertia of the wall with respect to
                              point O)
    :param alpha: angle between the vertical edge of the wall, passing through O, and the radius R 
    :param time_acceleration_interpolation: interpolation object time-accelarations
    :param gravity: gravity acceleration, i.e 9.81 m/s^2
    :return:
    """
    theta, omega = y

    #compute u_theta 
    #u_theta=0
    
    u_theta, p_theta_quadro = calculate_u_theta(theta)
    
    derivs = [omega, -p_theta_quadro * (np.sin(alpha - theta)-u_theta/R) + p_quadro * time_acceleration_interpolation(t) * np.cos(alpha - theta) / gravity]
    return derivs

def get_number_of_consecutive_negative_thetas(derivative_solution):
    negative_theta_values = 0
    if derivative_solution is not None:

        # iterate over the negative theta values and stop when a positive is found
        for i in range(len(derivative_solution[0][:, 0])):
            if derivative_solution[0][i, 0] <= 0:
                negative_theta_values += 1
            else:
                negative_theta_values -= 1
                break
    return negative_theta_values


def find_rotations(time_array, accelerations, calculation_accuracy):
    """
    Functions to calculate the rotations based on the time and the accelerations.

    :param time_array:
    :param accelerations:
    :param calculation_accuracy: gives the frequency of how many times the odeint function
        # should be called. If the accuracy is low, all negative theta values are skipped, and the computation is much faster.
        # if it is high, no negative theta values is skipped and the odeint function is called for every interval in time.
    :return:
    """
    rotations = []
    velocities = []
    skip_forward_list = []
    u_theta_real = []

    time_list_of_lists = [t]

    psoln = None

    activation = 0
    current_index = 0

    TEST_STOP = 1600

    while current_index < len(time_array) and current_index < TEST_STOP:
        velocity_after_impact = 0

        rotations_at_step = []
        velocities_at_step = []
        skip_forward_at_step = 0
        u_theta_real_step = []

        # the first time that we are here we consider velocity_impact_after == omega0
        if psoln is None:
            theta_after_impact = theta0
            velocity_after_impact = omega0
        else:
            # check the second rotation. If it is greater than zero, we consider all the
            # positive rotations found by the calculation of the derivatives
            theta_after_impact = 0
            if psoln[0][1, 0] > 0:
                activation = 1
                #u_theta=0
                # iterate over positive theta values and stop when a negative value is found
                for i in range(len(psoln[0][:, 0])):
                    if psoln[0][i, 0] >= 0:
                        current_index += 1

                        rotations_at_step.append(psoln[0][i, 0] * 180 / np.pi)
                        velocities_at_step.append(psoln[0][i, 1])
                        velocity_urto_meno = velocities_at_step[-1]

                        # the velocity after the impact is the last velocity multiplied by a coefficient
                        velocity_after_impact = r * velocity_urto_meno
                        u_theta_real_step.append(calculate_u_theta(psoln[0][i, 0])[0])

                    else:
                        break
                skip_forward_at_step=len(rotations_at_step)                   

            else:
                # check the new start index here
                activation = 0
                velocity_after_impact = 0

        # rotations at_step is empty if it is the first iteration (psol is None) or
        # the second rotation in the calculation of the derivatives is not positive
        rotations.append(rotations_at_step)
        velocities.append(velocities_at_step)
        u_theta_real.append(u_theta_real_step)
            
        y_curr = [theta_after_impact, velocity_after_impact]

        negative_theta_values = get_number_of_consecutive_negative_thetas(psoln)

        # we skip all the negative thetas
        if activation == 0:
            skip_forward=int(negative_theta_values * (1 - calculation_accuracy)) + 1
            current_index += skip_forward
            skip_forward_at_step=skip_forward

        next_time_array = time_array[current_index:]
        time_list_of_lists.append(next_time_array)
        
        skip_forward_list.append(skip_forward_at_step)

        #print(activation, current_index, next_time_array[0],skip_forward_at_step)

        # for interpolation we need at least two entries.
        if current_index >= len(time_array) - 1:
            return rotations, time_list_of_lists, skip_forward_list, u_theta_real
        
        else:
            interpol = interpolate_time_accelerations(next_time_array, accelerations[current_index:])
            psoln = odeint(calculate_derivatives, y_curr, next_time_array, args=(p_quadro, alpha, interpol, g), full_output=True)

    return rotations, time_list_of_lists, skip_forward_list, u_theta_real
    #return skip_forward_at_step

def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny+dy, maxy+dy)
    
def plot_rotations(time_list_of_lists, rotations_list_of_lists):
    """
    Plot the rotations.

    :param time_list_of_lists: list of list of interval of times
    :param rotations: list of list of rotations
    :return:
    """

    # assert len(time_list_of_lists) == len(rotations_list_of_lists)
    
    t_total = [0]
    rotation_total = [0]
    u_theta_total = [0]
    for ndv in range(len(rotations)):
        if not rotations[ndv]: # if rotations[ndv] is empty
            for i in range(0, skip_forward_list[ndv]):
                t_total+= [t_total[-1]+tInc]
                rotation_total+=[0]
                u_theta_total+=[0.3]
                
        else:
            t_total += list(ts[ndv][: len(rotations[ndv])])
            rotation_total += rotations[ndv]
            u_theta_total += u_theta_real[ndv]
    
    interpol = interpolate.interp1d(t, acc)
    
    # double axis plot   
    # as negative accelerations are of no interest, the plot focuses on positive accelerations only
    fig, ax1 = plt.subplots()
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Rotations $\Theta$ (deg)", color="red")
    plt1=ax1.plot(t_total, rotation_total, "r-", linewidth=2, label = 'Rotations')
    ax1.tick_params(axis="y", labelcolor="red")
    ax1.set_ylim([0, max(rotation_total) +0.3])
       
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.set_ylabel("Seismic acceleration [m/$s^2$]", color="blue")  # we already handled the x-label with ax1
    plt2=ax2.plot(t_total, interpol(t_total), linewidth = 0.6, label = 'Accelerations')
    plt3=ax2.plot(t_total, [acceleration_of_start_rocking] * len(t_total), "b--", linewidth = 1.5, label = 'Acceleration of activation')
    ax2.plot(t_total, [0] * len(t_total), "b", linewidth = 0.6)
    ax2.tick_params(axis="y", labelcolor="blue")
    #ax2.set_ylim([interpol(t_total).min(),interpol(t_total).max()+1])
    ax2.set_ylim([0,6])
    #add legend
    lns=plt1+plt2+plt3
    labs = [l.get_label() for l in lns]
    ax2.legend(lns, labs, loc=0)
    
    #progression of u_theta with theta
    fig, ax = plt.subplots()
    plt.plot
    ax.grid()
    #ax.plot(theta_array[theta_array>0],u_theta_array[u_theta_array<0.3])
    ax.plot(rotation_total[1:],u_theta_total[1:])
    ax.set(xlabel='$\Theta$ (deg)',ylabel='$U_{\Theta}$ (m)',title='complete dynamic of $\Theta$')
    
    plt.show()
    return t_total, rotation_total

b = 0.3  # mezza base
h = 2.4  # mezza altezza
H = 2 * h
zeta = h / b

# TODO: we consider only two dimensional wall
vol = 2 * b * 2 * h
density = 203
g = 9.81
m = vol * density
R = np.sqrt((b) ** 2 + (h) ** 2)
I = (4 / 3) * m * R ** 2
p_quadro = m * g * R / I
p = np.sqrt(p_quadro)
alpha = np.arctan(b / h)
acceleration_of_start_rocking=np.tan(alpha)*g #static acceleration of start rocking

# parameters for u_teta
B=2*b
l_h=1
k_n=6e6 #Pa
f_m=32739 #Pa 
teta_j0=2*m*g/(k_n*B**2*l_h)
a_c=2*m*g/(f_m*l_h)
teta_jc=f_m/(k_n*a_c)
list1 = []
list2 = []

# Make time array for solution
tStop = 8.000
tInc = 0.005
t = np.linspace(0.0, tStop, int(tStop / tInc + 1))


# Initial values
theta0 = 0 * np.pi / 180  # initial angular displacement RADIANTI
omega0 = 0.0  # initial angular velocity
#y_curr = [theta0, omega0]

e_1s = 1.05 * (1 - 2 * m * R ** 2 / I * np.sin(alpha) ** 2) ** 2 * (1 - 2 * m * R ** 2 / I * np.cos(alpha) ** 2)
r = e_1s

# Here everything starts

acc = read_accelerations("acceleration_laquila.xlsx",t)

interpol = interpolate_time_accelerations(t, acc)

rotations, ts, skip_forward_list, u_theta_real = find_rotations(t, acc, 0.7)

t_total, rotation_total = plot_rotations(ts, rotations)
