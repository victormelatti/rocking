# ground motion
# theta in function of time for accelerogram ground motion
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import xlrd
from scipy import interpolate
#from constants import HALF_BASE, HALF_HEIGHT, DEPTH, DENSITY, GRAVITY,vINTERFACE_STIFFNESS, COMPRESSIVE_STRENGTH, T_STOP, T_INC
import constants as cs

def compute_base_height(HALF_BASE, HALF_HEIGHT):
    BASE = HALF_BASE * 2
    HEIGHT = HALF_HEIGHT * 2
    return BASE, HEIGHT

def compute_volume(BASE, HEIGHT):
    volume = BASE * HEIGHT
    return volume

def compute_mass(VOLUME, DENSITY):
    mass = VOLUME * DENSITY
    return mass
    
def compute_radius (HALF_BASE, HALF_HEIGHT):
    radius = np.sqrt(HALF_BASE ** 2 + HALF_HEIGHT ** 2)    # meters
    return radius

def compute_polar_moment_inertia(MASS,radius):
    polar_moment_inertia = (4 / 3) * MASS * radius ** 2
    return polar_moment_inertia

def compute_alpha(HALF_BASE, HALF_HEIGHT):
    alpha = np.arctan(HALF_BASE / HALF_HEIGHT)
    return alpha

def compute_zeta(HALF_BASE, HALF_HEIGHT):
    zeta = HALF_HEIGHT / HALF_BASE
    return zeta

def compute_acceleration_of_start_rocking(alpha, GRAVITY):
    acceleration_of_start_rocking = np.tan(alpha) * GRAVITY  # static acceleration of start rocking
    return acceleration_of_start_rocking

def compute_velocity_reduction_coefficient (MASS, radius, alpha, I):
    velocity_reduction_coefficient= 1.05 * (1 - 2 * MASS * radius ** 2 / I * np.sin(alpha) ** 2) ** 2 * (1 - 2 * MASS * radius ** 2 / I * np.cos(alpha) ** 2)
    return velocity_reduction_coefficient

def compute_fully_compressed_theta(MASS, GRAVITY, INTERFACE_STIFFNESS, BASE, DEPTH):
    fully_compressed_theta = 2 * MASS * GRAVITY / (INTERFACE_STIFFNESS * BASE ** 2 * DEPTH)
    return fully_compressed_theta

def compute_fully_cracked_theta(COMPRESSIVE_STRENGTH, DEPTH, INTERFACE_STIFFNESS, MASS, GRAVITY):
    fully_cracked_theta = 0.5 * COMPRESSIVE_STRENGTH **2 * DEPTH / (INTERFACE_STIFFNESS * MASS * GRAVITY )
    return fully_cracked_theta

def compute_time_array(T_STOP, T_INC):
    time_array = np.linspace(0.0, T_STOP, int(T_STOP / T_INC + 1))
    return time_array

base, height = compute_base_height(cs.HALF_BASE, cs.HALF_HEIGHT)
volume = compute_volume(base, height)
mass = compute_mass(volume, cs.DENSITY)
radius = compute_radius(cs.HALF_BASE, cs.HALF_HEIGHT)
polar_moment_inertia = compute_polar_moment_inertia(mass,radius)
alpha = compute_alpha(cs.HALF_BASE, cs.HALF_HEIGHT)
zeta = compute_zeta(cs.HALF_BASE, cs.HALF_HEIGHT)

print(cs.INTERFACE_STIFFNESS)
acceleration_of_start_rocking = compute_acceleration_of_start_rocking(alpha, cs.GRAVITY)
velocity_reduction_coefficient = compute_velocity_reduction_coefficient(cs.MASS, radius, alpha, polar_moment_inertia)
fully_compressed_theta = compute_fully_compressed_theta(cs.MASS, cs.GRAVITY, cs.INTERFACE_STIFFNESS, base, cs.DEPTH)
fully_cracked_theta = compute_fully_cracked_theta(cs.COMPRESSIVE_STRENGTH, cs.DEPTH, cs.INTERFACE_STIFFNESS, mass, cs.GRAVITY)
time_array = compute_time_array(cs.T_STOP, cs.T_INC)


def read_accelerations(filename, time_array):
    """
    Function that reads accelerations from a source file

    :param filename: the name of source file
    :return: a list of accelerations
    """
    wb = xlrd.open_workbook(filename)
    s = wb.sheet_by_name("acceleration strong")
    time_array_length = len(time_array)
    accelerations = [s.cell_value(i, 4) * 0.002 for i in range(time_array_length)]

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


def calculate_u_theta(rotation):
    """
    Function that computes u_theta given a rotation

    :param rotation: rotation in radians
    :return: new utheta
    """

    u_theta = 0
    if rotation <= teta_j0:
        u_theta = B / 2 - (B ** 3 * k_n * l_h * rotation / (12 * m * g))

    if teta_j0 < rotation <= teta_jc:
        u_theta = (1 / 3) * np.sqrt(2 * m * g / (k_n * l_h * rotation))

    if rotation > teta_jc:
        u_theta = 0.5 * (m * g / (f_m * l_h) + (f_m ** 3 * l_h) / (12 * m * g * k_n ** 2 * rotation ** 2))

    return u_theta


def calculate_p_quadro(u_theta):
    """
    Function that computes p_quadro given u_theta
    :param u_theta:
    :return: p_quadro
    """

    r_theta = np.sqrt((R * np.cos(alpha)) ** 2 + (R * np.sin(alpha) - u_theta) ** 2)
    i_theta = m / 3 * R ** 2 + m * r_theta ** 2
    p_quadro = m * g * R / i_theta
    return p_quadro


def calculate_derivatives(y, t, R, alpha, time_acceleration_interpolation, gravity):
    """
    Function that scipy.odeint takes in input in order to calculate the
    rotation theta and the angular velocity omega.

    :param y: output of the differential equation
    :param t: input of the differential equation: array of equal-spaced float
    :param R: distance between the centroid G and the rotation point O
    :param alpha: angle between the vertical edge of the wall, passing through O, and the radius R 
    :param time_acceleration_interpolation: interpolation object time-accelarations
    :param gravity: gravity acceleration, i.e 9.81 m/s^2
    :return:
    """
    theta, omega = y

    u_theta = calculate_u_theta(theta)
    p_quadro = calculate_p_quadro(u_theta)

    derivs = [
        omega,
        -p_quadro * (np.sin(alpha - theta) - u_theta / R) + p_quadro * time_acceleration_interpolation(t) * np.cos(alpha - theta) / gravity,
    ]
    return derivs


def get_number_of_consecutive_non_positive_thetas(derivative_solution):
    negative_theta_values = 0
    if derivative_solution is not None:

        # iterate over the non-positive theta values and stop when a positive is found
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

                # iterate over positive theta values and stop when a negative value is found
                for i in range(len(psoln[0][:, 0])):
                    if psoln[0][i, 0] >= 0:
                        current_index += 1

                        rotations_at_step.append(psoln[0][i, 0] * 180 / np.pi)
                        velocities_at_step.append(psoln[0][i, 1])
                        velocity_urto_meno = velocities_at_step[-1]

                        # the velocity after the impact is the last velocity multiplied by a coefficient
                        velocity_after_impact = e_1s * velocity_urto_meno
                        u_theta_real_step.append(calculate_u_theta(psoln[0][i, 0]))

                    else:
                        break
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

        non_positive_theta_values = get_number_of_consecutive_non_positive_thetas(psoln)

        # we skip all the non-positive thetas
        if activation == 0:
            skip_forward = int(non_positive_theta_values * (1 - calculation_accuracy)) + 1
            current_index += skip_forward

        next_time_array = time_array[current_index:]
        time_list_of_lists.append(next_time_array)

        # for interpolation we need at least two entries.
        if current_index >= len(time_array) - 1:
            return rotations, time_list_of_lists, u_theta_real

        else:
            interpol = interpolate_time_accelerations(next_time_array, accelerations[current_index:])
            psoln = odeint(calculate_derivatives, y_curr, next_time_array, args=(R, alpha, interpol, g), full_output=True)

    return rotations, time_list_of_lists, u_theta_real


def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1 - y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny + dy, maxy + dy)


def line_chart(x, y, color, label_x, label_y):
    # TODO: use this function
    pass


def plot_rotations(time_list_of_lists, rotations_list_of_lists, accelerations):
    """
    Plot the rotations.

    :param time_list_of_lists: list of list of interval of times
    :param rotations: list of list of rotations
    :return:
    """

    t_total = [0]
    rotation_total = [0]
    u_theta_total = [0]

    for ndv in range(len(rotations)):
        if not rotations[ndv]:  # if rotations[ndv] is empty
            t_total.append(t_total[-1] + tInc)
            rotation_total.append(0)
            u_theta_total.append(b)
        else:
            t_total += list(ts[ndv][: len(rotations[ndv])])
            rotation_total += rotations[ndv]
            u_theta_total += u_theta_real[ndv]

    interpol = interpolate_time_accelerations(t, accelerations)

    # double axis plot
    # as negative accelerations are of no interest, the plot focuses on positive accelerations only
    fig, ax1 = plt.subplots()
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Rotations $\Theta$ (deg)", color="red")
    plt1 = ax1.plot(t_total, rotation_total, "r-", linewidth=2, label="Rotations")
    ax1.tick_params(axis="y", labelcolor="red")
    ax1.set_ylim([0, max(rotation_total) + 0.3])

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.set_ylabel("Seismic acceleration [m/$s^2$]", color="blue")  # we already handled the x-label with ax1
    plt2 = ax2.plot(t_total, interpol(t_total), linewidth=0.6, label="Accelerations")
    plt3 = ax2.plot(t_total, [acceleration_of_start_rocking] * len(t_total), "b--", linewidth=1.5, label="Acceleration of activation")
    ax2.plot(t_total, [0] * len(t_total), "b", linewidth=0.6)
    ax2.tick_params(axis="y", labelcolor="blue")
    # ax2.set_ylim([interpol(t_total).min(),interpol(t_total).max()+1])
    ax2.set_ylim([0, 6])
    # add legend
    lns = plt1 + plt2 + plt3
    labs = [l.get_label() for l in lns]
    ax2.legend(lns, labs, loc=0)

    # progression of u_theta with theta
    fig, ax = plt.subplots()
    ax.grid()
    # ax.plot(theta_array[theta_array>0],u_theta_array[u_theta_array<0.3])
    ax.plot(rotation_total[1:], u_theta_total[1:])
    ax.set(xlabel="$\Theta$ (deg)", ylabel="$U_{\Theta}$ (m)", title="complete dynamic of $\Theta$")

    plt.show()
    return t_total, rotation_total

# Initial values
theta0 = 0 * np.pi / 180  # initial angular displacement in radians
omega0 = 0.0 # initial angular velocity

# Here everything starts

acc = read_accelerations("acceleration_laquila.xlsx", t)

interpol = interpolate_time_accelerations(t, acc)

rotations, ts, u_theta_real = find_rotations(t, acc, 1)

t_total, rotation_total = plot_rotations(ts, rotations, acc)
