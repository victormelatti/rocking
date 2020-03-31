import numpy as np
import matplotlib.pyplot as plt
import xlrd
from scipy import interpolate
from scipy.integrate import odeint
import constants_anchor as cs


def compute_base_height(half_base, half_height):
    base = half_base * 2
    height = half_height * 2
    return base, height


def compute_volume(base, height):
    volume = base * height
    return volume


def compute_mass(volume, density):
    mass = volume * density
    return mass


def compute_radius(half_base, half_height):
    radius = np.sqrt(half_base ** 2 + half_height ** 2)  # meters
    return radius


def compute_polar_moment_inertia(mass, radius):
    polar_moment_inertia = (4 / 3) * mass * radius ** 2
    return polar_moment_inertia


def compute_alpha(half_base, half_height):
    alpha = np.arctan(half_base / half_height)
    return alpha


def compute_zeta(half_base, half_height):
    zeta = half_height / half_base
    return zeta


def compute_acceleration_of_start_rocking(alpha, gravity):
    acceleration_of_start_rocking = np.tan(alpha) * gravity  # static acceleration of start rocking
    return acceleration_of_start_rocking


def compute_velocity_reduction_coefficient(mass, radius, alpha, I):
    velocity_reduction_coefficient = (
        1.05 * (1 - 2 * mass * radius ** 2 / I * np.sin(alpha) ** 2) ** 2 * (1 - 2 * mass * radius ** 2 / I * np.cos(alpha) ** 2)
    )
    return velocity_reduction_coefficient


def compute_fully_compressed_theta(mass, gravity, interface_stiffness, base, depth):
    theta_fully_compressed = 2 * mass * gravity / (interface_stiffness * base ** 2 * depth)
    return theta_fully_compressed


def compute_fully_cracked_theta(compressive_strength, depth, interface_stiffness, mass, gravity):
    theta_fully_cracked = 0.5 * compressive_strength ** 2 * depth / (interface_stiffness * mass * gravity)
    return theta_fully_cracked


def compute_time_array(t_stop, t_inc):
    time_array = np.linspace(0.0, t_stop, int(t_stop / t_inc + 1))
    return time_array

#anchor functions

def compute_anchor_height(height, distance_from_top):
    anchor_height = height - distance_from_top
    return anchor_height

def compute_young_modulus(sigma_yielding, epsilon_yielding):
    young_modulus = sigma_yielding / epsilon_yielding
    return young_modulus

def compute_anchor_area(number_of_anchors, anchor_diameter):
    anchor_area = number_of_anchors * np.pi * anchor_diameter ** 2 / 4 #m^2
    return anchor_area

def compute_anchor_yielding_force(sigma_yielding, anchor_area):
    anchor_yielding_force = sigma_yielding * anchor_area
    return anchor_yielding_force

def rotation_of_yielding(epsilon_yielding, anchor_length, anchor_height):
    rotation_yielding = np.arctan(epsilon_yielding*anchor_length/anchor_height)
    return rotation_yielding

def rotation_of_ultimate_deformation(epsilon_ultimate, anchor_length, anchor_height):
    rotation_ultimate_deformation = np.arctan(epsilon_ultimate * anchor_length/anchor_height)
    return rotation_ultimate_deformation

def calculate_rotation_max_so_far(rotation, rotation_max):  
    """
    computes maximum rotation reached so far in the rocking motion of the block
    """
    if rotation >= rotation_max:
        rotation_max = rotation        
    return rotation_max

def calculate_anchor_force(
        rotation,
        rotation_max,  
        rotation_residual,
        rotation_yielding,
        rotation_ultimate_deformation,
        anchor_yielding_force,
        ):
    """"
    function that computes the force in a horizontal steel tie placed at anchor_height from the ground.
    
    """

    rotation_max = calculate_rotation_max_so_far(rotation, rotation_max)  
    
    anchor_force_max_so_far = anchor_yielding_force* (rotation_max - rotation_residual)/rotation_yielding

    #zero branch
    if  rotation <= rotation_residual:
        anchor_force = 0   
    
    #loading and uloading branch
    if rotation_residual < rotation <= rotation_max: 
        anchor_force = anchor_force_max_so_far * (1 + (rotation - rotation_max) / (rotation_max - rotation_residual))

    #horizontal branch
    if rotation_yielding + rotation_residual <= rotation <= rotation_ultimate_deformation: 
        anchor_force = anchor_yielding_force        
        rotation_residual = rotation_max - rotation_yielding

    #rupture branch
    if rotation > rotation_ultimate_deformation:
        anchor_force = 0 
        rotation_residual = rotation_max

    
    return anchor_force, rotation_max, rotation_residual


def read_accelerations(filename, time_array):
    """
    Function that reads accelerations from a source file

    :param filename: the name of source file
    :return: a list of accelerations
    """
    open_excel_file = xlrd.open_workbook(filename)
    excel_sheet = open_excel_file.sheet_by_name("acceleration strong")
    time_array_length = len(time_array)
    accelerations = [excel_sheet.cell_value(i, 4) * 0.01 for i in range(time_array_length)]

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


def calculate_rotation_point(
    rotation, theta_fully_compressed, theta_fully_cracked, base, interface_stiffness, compressive_strength, depth, gravity, mass
):
    """
    Function that computes u_theta given a rotation

    :param rotation: rotation in radians
    :return: new utheta
    """

    rotation_point = 0
    if rotation <= theta_fully_compressed:
        rotation_point = base / 2 - (base ** 3 * interface_stiffness * depth * rotation / (12 * mass * gravity))

    if theta_fully_compressed < rotation <= theta_fully_cracked:
        rotation_point = (1 / 3) * np.sqrt(2 * mass * gravity / (interface_stiffness * depth * rotation))

    if rotation > theta_fully_cracked:
        rotation_point = 0.5 * (
            mass * gravity / (compressive_strength * depth)
            + (compressive_strength ** 3 * depth) / (12 * mass * gravity * interface_stiffness ** 2 * rotation ** 2)
        )

    return rotation_point


def calculate_frequency_parameter(rotation_point, radius, mass, gravity, alpha):
    """
    Function that computes p_quadro given u_theta
    :param u_theta:
    :return: p_quadro
    """

    radius_theta = np.sqrt((radius * np.cos(alpha)) ** 2 + (radius * np.sin(alpha) - rotation_point) ** 2)
    polar_moment_inertia_theta = mass / 3 * radius_theta ** 2 + mass * radius_theta ** 2
    frequency_parameter = mass * gravity * radius / polar_moment_inertia_theta

    return frequency_parameter

rotation_max, rotation_residual = (0, 0)

def calculate_derivatives(
    y,
    time_array,
    radius,
    alpha,
    time_acceleration_interpolation,
    gravity,
    theta_fully_compressed,
    theta_fully_cracked,
    base,
    mass,
    interface_stiffness,
    compressive_strength,
    depth,
    anchor_height,
    rotation_max,
    rotation_residual,
    rotation_yielding, 
    rotation_ultimate_deformation,
    anchor_yielding_force
):
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

    rotation_point = calculate_rotation_point(
        theta, theta_fully_compressed, theta_fully_cracked, base, interface_stiffness, compressive_strength, depth, gravity, mass
    )

    frequency_parameter = calculate_frequency_parameter(rotation_point, radius, mass, gravity, alpha)

    anchor_force, rotation_max, rotation_residual = calculate_anchor_force( theta, rotation_max, rotation_residual, rotation_yielding, rotation_ultimate_deformation, anchor_yielding_force)
    
    derivs = [
        omega,
        - frequency_parameter * (np.sin(alpha - theta) - rotation_point / radius)
        + frequency_parameter * time_acceleration_interpolation(time_array) * np.cos(alpha - theta) / gravity 
        - frequency_parameter * anchor_force * anchor_height * np.cos(theta) / (mass * radius * gravity),
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



def find_rotations(
    time_array,
    accelerations,
    calculation_accuracy,
    theta_fully_compressed,
    theta_fully_cracked,
    base,
    interface_stiffness,
    compressive_strength,
    depth,
    gravity,
    mass,
    radius,
    alpha,
    theta0,
    omega0,
    velocity_reduction_coefficient,
    anchor_height,
    rotation_yielding, 
    rotation_ultimate_deformation,
    anchor_yielding_force
):
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
    rotation_points = []
    anchor_force = []
    rotation_max, rotation_residual = (0, 0)

    time_list_of_lists = [time_array]

    psoln = None
    current_index = 0

    TEST_STOP = 1600
    
    rotation_point_i = base / 2

    while current_index < len(time_array) and current_index < TEST_STOP:

        theta_after_impact = theta0 if psoln is None else 0
        velocity_after_impact = omega0 if psoln is None else 0
        activation = 0

        rotations_at_step = []
        rotation_points_step = []
        anchor_force_at_step = []

        if psoln is not None and psoln[0][1, 0] > 0:
            activation = 1
            i = 0

            # iterate over positive theta values and stop when a negative value is found
            while i < len(psoln[0][:, 0]) and psoln[0][i, 0] >= 0:
                current_index += 1
                rotation = psoln[0][i, 0] * 180 / np.pi
                velocity = psoln[0][i, 1]

                rotations_at_step.append(rotation)
                velocity_after_impact = velocity_reduction_coefficient * velocity

                rotation_point_i = calculate_rotation_point(
                    psoln[0][i, 0],
                    theta_fully_compressed,
                    theta_fully_cracked,
                    base,
                    cs.INTERFACE_STIFFNESS,
                    cs.COMPRESSIVE_STRENGTH,
                    cs.DEPTH,
                    cs.GRAVITY,
                    mass,
                )
                rotation_points_step.append(rotation_point_i)
                                                
                anchor_force_i, rotation_max, rotation_residual = calculate_anchor_force(
                    psoln[0][i, 0],
                    rotation_max,
                    rotation_residual,
                    rotation_yielding,
                    rotation_ultimate_deformation,
                    anchor_yielding_force)
            
                anchor_force_at_step.append(anchor_force_i)
                i += 1
        
        rotations.append(rotations_at_step)
        rotation_points.append(rotation_points_step)
        anchor_force.append(anchor_force_at_step)

        y_curr = [theta_after_impact, velocity_after_impact]

        non_positive_theta_values = get_number_of_consecutive_non_positive_thetas(psoln)

        # we skip all the non-positive thetas
        if activation == 0:
            skip_forward = int(non_positive_theta_values * (1 - calculation_accuracy)) + 1
            current_index += skip_forward

        next_time_array = time_array[current_index:]
        time_list_of_lists.append(next_time_array)
        print(current_index, activation)    

        # for interpolation we need at least two entries.
        if current_index >= len(time_array) - 1:
            return rotations, time_list_of_lists, rotation_points, anchor_force

        else:

            interpol = interpolate_time_accelerations(next_time_array, accelerations[current_index:])
            psoln = odeint(
                calculate_derivatives,
                y_curr,
                next_time_array,
                args=(
                    radius,
                    alpha,
                    interpol,
                    cs.GRAVITY,
                    theta_fully_compressed,
                    theta_fully_cracked,
                    base,
                    mass,
                    cs.INTERFACE_STIFFNESS,
                    cs.COMPRESSIVE_STRENGTH,
                    cs.DEPTH,
                    anchor_height,
                    rotation_max,
                    rotation_residual,
                    rotation_yielding, 
                    rotation_ultimate_deformation,
                    anchor_yielding_force
                ),
                full_output=True,
            )
            
    return rotations, time_list_of_lists, rotation_points, anchor_force


def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1 - y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny + dy, maxy + dy)


def line_chart(x, y, color, label_x, label_y):
    # TODO: use this function TODO: to do what?
    pass


def plot_rotations(
    time_list_of_lists,
    rotations_list_of_lists,
    rotation_point_list_of_lists,
    anchor_force_list_of_lists,
    accelerations,
    time_increment,
    half_base,
    time_array,
    acceleration_of_start_rocking,
):
    """
    Plot the rotations.

    :param time_list_of_lists: list of list of interval of times
    :param rotations: list of list of rotations
    :return:
    """

    t_total = [0]
    rotation_total = [0]
    rotation_point_total = [0]
    anchor_force_total = [0]

    for ndv in range(len(rotations_list_of_lists)):
        if not rotations_list_of_lists[ndv]:  # if rotations_list_of_lists[ndv] is empty
            t_total.append(t_total[-1] + time_increment)
            rotation_total.append(0)
            rotation_point_total.append(half_base)
            anchor_force_total.append(0)
        else:
            t_total += list(time_list_of_lists[ndv][: len(rotations_list_of_lists[ndv])])
            rotation_total += rotations_list_of_lists[ndv]
            rotation_point_total += rotation_point_list_of_lists[ndv]
            anchor_force_total += anchor_force_list_of_lists[ndv]

    interpol = interpolate_time_accelerations(time_array, accelerations)

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
    ax2.set_ylim([0, max(interpol(t_total)) + 0.5])
    # add legend
    lns = plt1 + plt2 + plt3
    labs = [l.get_label() for l in lns]
    ax2.legend(lns, labs, loc=0)

    # progression of u_theta with theta
    fig, ax = plt.subplots()
    ax.grid()
    ax.plot(rotation_total[1:], rotation_point_total[1:])
    ax.set(xlabel="$\Theta$ (deg)", ylabel="$U_{\Theta}$ (m)", title="$U_{\Theta}$ vs $\Theta$")

    # progression of anchor_force with theta
    fig, ax = plt.subplots()
    ax.grid()
    ax.plot(rotation_total, anchor_force_total)
    ax.set(xlabel="$\Theta$ (deg)", ylabel="Anchor force (KN)", title=" $\Theta$")


    plt.show()
    return t_total, rotation_total, anchor_force_total     
        