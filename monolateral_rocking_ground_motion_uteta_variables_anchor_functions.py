import constants_anchor as cs
import functions_anchor as fn


__name__ = "__main__"
if __name__ == "__main__":
    print("This file solves the equation of motion of a rocking block with anchor")

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
    
    anchor_height = fn.compute_anchor_height(height, cs.DISTANCE_FROM_TOP)
    
    young_modulus = fn.compute_young_modulus(cs.SIGMA_YIELDING, cs.EPSILON_YIELDING)

    anchor_area = fn.compute_anchor_area(cs.NUMBER_OF_ANCHORS, cs.ANCHOR_DIAMETER)
    
    anchor_yielding_force = fn.compute_anchor_yielding_force(cs.SIGMA_YIELDING, anchor_area)
    
    rotation_yielding = fn.rotation_of_yielding(cs.EPSILON_YIELDING, cs.ANCHOR_LENGTH, anchor_height)
    
    rotation_ultimate_deformation = fn.rotation_of_ultimate_deformation(cs.EPSILON_ULTIMATE, cs.ANCHOR_LENGTH, anchor_height)

    rotations, ts, rotation_point, anchor_force = fn.find_rotations(
        time_array,
        accelerations,
        0.5,
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
        anchor_height,
        rotation_yielding, 
        rotation_ultimate_deformation,
        anchor_yielding_force
    )

    t_total, rotation_total, anchor_force_total = fn.plot_rotations(
        ts, rotations, rotation_point, anchor_force, accelerations, cs.T_INC, cs.HALF_BASE, time_array, acceleration_of_start_rocking
    )

    #m_anchor = fn.plot_moments(moments)

