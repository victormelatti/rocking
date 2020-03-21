# ground motion
# theta in function of time for accelerogram ground motion
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import xlrd
from scipy import interpolate


def read_accelerations(filename):
    """
    Function that reads accelerations from a source file

    :param filename: the name of source file
    :return: a list of accelerations
    """
    # seismic action
    wb = xlrd.open_workbook(filename)
    s = wb.sheet_by_name("acceleration strong")

    # time=np.linspace(0., tStop, int(tStop/tInc)+1)
    # time=time[6000:10000]
    accelerations = [s.cell_value(i, 4) * 0.01 for i in range(s.nrows)]

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


def calculate_derivatives(y, t, p_quadro, alpha, time_acceleration_interpolation, gravity):
    """
    Function that scipy.odeint takes in input in order to calculate the
    rotation theta and the angular velocity omega.

    :param y: output of the differential equation
    :param t: input of the differential equation: array of equal-spaced float
    :param p_quadro: TODO: write here what p_quadro is
    :param alpha: TODO: write what alpha is
    :param time_acceleration_interpolation: interpolation object time-accelarations
    :param gravity: gravity acceleration, i.e 9.81 m/s^2
    :return:
    """
    theta, omega = y
    # p_quadro, alpha, interpol, g = params

    derivs = [omega, -p_quadro * (np.sin(alpha - theta)) + p_quadro * time_acceleration_interpolation(t) * np.cos(alpha - theta) / gravity]
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


def find_rotations(time_array, accelerations, calculation_accuracy=0.7):
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

    time_list_of_lists = [t]

    psoln = None

    activation = 0
    current_index = 0

    TEST_STOP = 1600

    while current_index < len(time_array) and current_index < TEST_STOP:
        velocity_after_impact = 0

        rotations_at_step = []
        velocities_at_step = []

        # the first time that we are here we consider velocity_impact_after == omega0
        if psoln is None:
            velocity_after_impact = omega0
        else:
            # check the second rotation. If it is greater than zero, we consider all the
            # positive rotations found by the calculation of the derivatives
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
                        velocity_after_impact = r * velocity_urto_meno
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

        y_curr = [0, velocity_after_impact]

        negative_theta_values = get_number_of_consecutive_negative_thetas(psoln)

        # we skip all the negative thetas
        if activation == 0:
            current_index += int(negative_theta_values * (1 - calculation_accuracy)) + 1

        print(activation, current_index)

        next_time_array = time_array[current_index:]
        time_list_of_lists.append(next_time_array)

        # for interpolation we need at least two entries.
        if current_index >= len(time_array) - 1:
            return rotations, time_list_of_lists
        else:
            interpol = interpolate_time_accelerations(next_time_array, accelerations[current_index:])
            psoln = odeint(calculate_derivatives, y_curr, next_time_array, args=(p_quadro, alpha, interpol, g), full_output=True)

    return rotations, time_list_of_lists


def plot_rotations(time_list_of_lists, rotations_list_of_lists):
    """
    Plot the rotations.

    :param time_list_of_lists: list of list of interval of times
    :param rotations: list of list of rotations
    :return:
    """

    # assert len(time_list_of_lists) == len(rotations_list_of_lists)

    # plot all together
    t_total = []
    rotation_total = []
    for ndv in range(len(rotations_list_of_lists)):
        t_total += list(time_list_of_lists[ndv][: len(rotations_list_of_lists[ndv])])
        rotation_total += rotations_list_of_lists[ndv]

    initial_part_time = list(np.linspace(0.0, t_total[0], int(t_total[0] / tInc + 1)))
    initial_part_rotations = [0] * len(initial_part_time)

    t_total = initial_part_time + t_total
    rotation_total = initial_part_rotations + rotation_total
    interpol = interpolate.interp1d(t, acc)

    fig, ax = plt.subplots()
    axes = plt.gca()
    ax.grid()
    ax.plot(t_total, rotation_total, "r-", linewidth=3)
    ax.plot(t_total, interpol(t_total))
    ax.set(xlabel="time (sec)", ylabel="$\Theta$ (deg)", title="complete dynamic of $\Theta$")
    axes = plt.gca()
    axes.set_xlim([0, 4])
    # axes.set_ylim([-1,10])
    # fig.savefig("Plot theta as a function of time_WITH IMPACT DISSIPATION_igor_0.5.png", dpi = 500, bbox_inches='tight')

    # plot su diversi axes
    # con questo plot si vede meglio come le accelerezioni positive producono rotazioni positive
    fig, ax1 = plt.subplots()

    ax1.set_xlabel("time (s)")
    ax1.set_ylabel("rotations $\Theta$ (deg)", color="red")
    ax1.plot(t_total, rotation_total, "r-", linewidth=3)
    ax1.tick_params(axis="y", labelcolor="red")

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.set_ylabel("seismic acceleration", color="blue")  # we already handled the x-label with ax1
    ax2.plot(t_total, interpol(t_total))
    ax2.plot(t_total, [0] * len(t_total), "b")
    ax2.tick_params(axis="y", labelcolor="blue")

    plt.show()


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


# Make time array for solution
tStop = 8.000
tInc = 0.005
t = np.linspace(0.0, tStop, int(tStop / tInc + 1))


# Initial values
theta0 = 0 * np.pi / 180  # initial angular displacement RADIANTI
omega0 = 0.0  # initial angular velocity
y_curr = [theta0, omega0]

e_1s = 1.05 * (1 - 2 * m * R ** 2 / I * np.sin(alpha) ** 2) ** 2 * (1 - 2 * m * R ** 2 / I * np.cos(alpha) ** 2)
r = e_1s


# Here everything starts

acc = read_accelerations("acceleration_laquila.xlsx")

interpol = interpolate_time_accelerations(t, acc)

rotations, ts = find_rotations(t, acc)

plot_rotations(ts, rotations)
