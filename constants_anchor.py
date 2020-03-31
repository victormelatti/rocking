import numpy as np
#geometrical parameters of the wall:


HALF_BASE = 0.3      # meters

HALF_HEIGHT = 2.4    # meters

DEPTH = 1.0   # meters

#physical parameters of the wall:
DENSITY = 203.0    # Kg/m^3

INTERFACE_STIFFNESS = 6e6  # Pa

COMPRESSIVE_STRENGTH = 32739  # Pa


T_STOP = 8.000   #sec

T_INC = 0.005    #sec

GRAVITY = 9.81    #m/s^2

THETA0 = 0 * np.pi / 180  # initial angular displacement in radians
OMEGA0 = 0.0 # initial angular velocity

#parameters for the anchor
ANCHOR_LENGTH = 4    #m
DISTANCE_FROM_TOP = 0.3
EPSILON_YIELDING = 0.002
EPSILON_ULTIMATE = 0.016
SIGMA_YIELDING = 235e3    #KN/m^2

#NUMBER_OF_ANCHORS = 9
NUMBER_OF_ANCHORS = float(input('Enter Number of anchors: '))
ANCHOR_DIAMETER = 0.012    #m


