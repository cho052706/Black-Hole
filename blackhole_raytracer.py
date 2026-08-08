import numpy as np

RS = 1.0
INCLINATION_DEG = 45.0
CAMERA_DIST = 16.0 * RS

WIDTH = 640
HEIGHT = 480

FOV_DEG = 42.0

incl = np.radians(INCLINATION_DEG)
cam_pos = np.array([np.sin(incl), 0.0, np.cos(incl)]) * CAMERA_DIST

forward = -cam_pos / np.linalg.norm(cam_pos)
world_up = np.array([0.0, 0.0, 1.0])
right = np.cross(forward, world_up)
right /= np.linalg.norm(right)
up = np.cross(right, forward)

fov = np.radians(FOV_DEG)
aspect = WIDTH / HEIGHT
xs = np.linspace(-1, 1, WIDTH) * np.tan(fov/2) * aspect
ys = np.linspace(1, -1, HEIGHT) * np.tan(fov/2)
xx, yy = np.meshgrid(xs, ys)

dirs = (
    forward[None, None, :]
    + xx[..., None] * right[None, None, :]
    + yy[..., None] * up[None, None, :]
)
dirs = dirs / np.linalg.norm(dirs, axis=-1, keepdims=True)
dirs = dirs.reshape(-1, 3)
N = dirs.shape[0]
pos0 = np.title(cam_pos, (N, 1))