import constants as cs
import functions as fn




if __name__ == "__main__":
    print("Hello")

    time_array = fn.compute_time_array(cs.T_STOP, cs.T_INC)

    accelerations = fn.read_accelerations("acceleration_laquila.xlsx", time_array)

    base, height = fn.compute_base_height(cs.HALF_BASE, cs.HALF_HEIGHT)

    volume = fn.compute_volume(base, height)

    mass = fn.compute_mass(volume, cs.DENSITY)

    radius = fn.compute_radius(cs.HALF_BASE, cs.HALF_HEIGHT)

    polar_moment_inertia = fn.compute_polar_moment_inertia(mass, radius)

    alpha = fn.compute_alpha(cs.HALF_BASE, cs.HALF_HEIGHT)

    zeta = fn.compute_zeta(cs.HALF_BASE, cs.HALF_HEIGHT)

    acceleration_of_start_rocking = fn.compute_acceleration_of_start_rocking(alpha, cs.GRAVITY)

    velocity_reduction_coefficient = fn.compute_velocity_reduction_coefficient(mass, radius, alpha, polar_moment_inertia)

    theta_fully_compressed = fn.compute_fully_compressed_theta(mass, cs.GRAVITY, cs.INTERFACE_STIFFNESS, base, cs.DEPTH)

    theta_fully_cracked = fn.compute_fully_cracked_theta(cs.COMPRESSIVE_STRENGTH, cs.DEPTH, cs.INTERFACE_STIFFNESS, mass, cs.GRAVITY)

    interpol = fn.interpolate_time_accelerations(time_array, accelerations)

    rotations, ts, rotation_point = fn.find_rotations(
        time_array,
        accelerations,
        1,
        theta_fully_compressed,
        theta_fully_cracked,
        base,
        cs.INTERFACE_STIFFNESS,
        cs.COMPRESSIVE_STRENGTH,
        cs.DEPTH,
        cs.GRAVITY,
        mass,
        radius,
        alpha,
        cs.THETA0,
        cs.OMEGA0,
        velocity_reduction_coefficient,
    )

    t_total, rotation_total = fn.plot_rotations(
        ts, rotations, rotation_point, accelerations, cs.T_INC, cs.HALF_BASE, time_array, acceleration_of_start_rocking
    )
