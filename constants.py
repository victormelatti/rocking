import numpy as np

HALF_BASE = 0.3      # meters
HALF_HEIGHT = 2.4    # meters
HEIGHT = 2 * HALF_HEIGHT
BASE = 2 * HALF_BASE

VOLUME = BASE * HEIGHT
DENSITY = 203
MASS = VOLUME * DENSITY

R = np.sqrt(HALF_BASE ** 2 + HALF_HEIGHT ** 2)
I = (4 / 3) * MASS * R ** 2
alpha = np.arctan(HALF_BASE / HALF_HEIGHT)