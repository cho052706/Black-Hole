import numpy as np

# Constants
RS = 1.0
INCLINATION_DEG = 45.0
CAMERA_DIST = 16.0 * RS

WIDTH = 640
HEIGHT = 480

FOV_DEG = 42.0
MAX_WINDINGS = 3.0
DPHI = 0.01

# Camera setup
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
pos0 = np.tile(cam_pos, (N, 1))

# Creating the orbital plane for the photon
r0 = np.linalg.norm(pos0, axis=1)
e_r = pos0 / r0[:, None]

normal = np.cross(e_r, dirs)
normlen = np.maximum(np.linalg.norm(normal, axis=1), 1e-12)
normal_hat = normal / normlen[:, None]
e_phi = np.cross(normal_hat, e_r)

cospsi = -np.einsum('ij,ij->i', dirs, e_r)
sinpsi = np.einsum('ij,ij->i', dirs, e_phi)
sinpsi = np.clip(sinpsi, 1e-6, 1.0)

u0 = 1.0 / r0
metric0 = np.clip(1.0 - RS * u0, 1e-12, None)
b = r0 * sinpsi / np.sqrt(metric0)
w0 = u0 * np.sqrt(metric0) * (cospsi / sinpsi)

u = u0.copy()
w = w0.copy()
phi = np.zeros(N)

active = np.ones(N, dtype=bool)
hit_disk = np.zeros(N, dtype=bool)

prev_z = pos0[:, 2].copy()
prev_pos = pos0.copy()


# RK4 integration
def deriv(u_, w_):
    return w_, -u_ + 1.5 * RS * (u_**2)

max_phi_total = MAX_WINDINGS * 2 * np.pi
nsteps = int(max_phi_total / DPHI)

for step in range(nsteps):
    if not active.any():
        break

    idx = active
    idxs = np.where(idx)[0]

    u_a, w_a = u[idx], w[idx]
    k1u, k1w = deriv(u_a, w_a)
    k2u, k2w = deriv(u_a + 0.5 * DPHI * k1u, w_a + 0.5 * DPHI * k1w)
    k3u, k3w = deriv(u_a + 0.5 * DPHI * k2u, w_a + 0.5 * DPHI * k2w)
    k4u, k4w = deriv(u_a + DPHI * k3u, w_a + DPHI * k3w)

    u_new = u_a + (DPHI / 6.0) * (k1u + 2*k2u + 2*k3u + k4u)
    w_new = w_a + (DPHI / 6.0) * (k1w + 2*k2w + 2*k3w + k4w)

    u[idx] = u_new
    w[idx] = w_new
    phi[idx] += DPHI

    r_new = 1.0 / np.clip(u_new, 1e-8, None)
    cosphi = np.cos(phi[idx])
    sinphi = np.sin(phi[idx])
    pos_new = r_new[:, None] * (cosphi[:, None] * e_r[idx] + sinphi[;, None] * e_phi[idx])
    z_new = pos_new[:, 2]