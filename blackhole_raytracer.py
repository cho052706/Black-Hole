import numpy as np

RS = 1.0
INCLINATION_DEG = 45.0
CAMERA_DIST = 16.0 * RS

incl = np.radians(INCLINATION_DEG)
cam_pos = np.array([np.sin(incl), 0.0, np.cos(incl)]) * CAMERA_DIST

forward = -cam_pos / np.linalg.norm(cam_pos)
world_up = np.array([0.0, 0.0, 1.0])
right = np.cross(forward, world_up)
right /= np.linalg.norm(right)
up = np.cross(right, forward)