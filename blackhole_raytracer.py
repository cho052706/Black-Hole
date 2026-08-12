import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter

# Constants
# To prevent a funny artifact use DPHI = 0.003, WIDTH = 960, HEIGHT = 720
RS = 1.0  # Schwarzschild radius
INCLINATION_DEG = 88.0  # Angle from the z-axis
TILT_DEG = -5.0  # Angle of camera rotation (clockwise is positive)
CAMERA_DIST = 30.0 * RS  # Distance of camera from the black hole

# Output resolution
WIDTH = 960
HEIGHT = 720

FOV_DEG = 42.0  # Camera field of view (degrees)
DPHI = 0.003  # RK4 step size in the swept angle phi
ESCAPE_R = 80.0 * RS  # Distance where a ray becomes "escaped"

# How many times a ray can loop around (impacts the clarity of photon sphere)
MAX_WINDINGS = 5.0

# Accretion disk dimensions
DISK_INNER = 3.0 * RS
DISK_OUTER = 16.0 * RS

BLOOM_ENABLED = True  # Adds a glow to specified pixels

DISK_BLOOM_THRESHOLD = 0.05  # Only pixels brighter than this glow
DISK_BLOOM_BOOST = 2.5  # Multiplies the light before blurring
# Formatted as (blur radius in pixels, contribution weight)
DISK_BLOOM_LAYERS = [(7.0, 0.9), (15.0, 0.5), (40.0, 0.3)]

STAR_BLOOM_THRESHOLD = 0.50
STAR_BLOOM_BOOST = 1.5
STAR_BLOOM_LAYERS = [(7.0, 0.9), (15.0, 0.5), (40.0, 0.3)]

# Star paramiters
STAR_DENSITY = 0.0015  # Probability that a bakground pixel will be a star
SEED = 7

rng = np.random.default_rng(SEED)


def deriv(u_, w_):
    """
    Derivatives for the first-order system equivalent to the Schwarzschild
    photon orbit equation d^2u/dphi^2 = -u + (3/2) RS u^2

    Parameters:
        u_: the current value(s) of u = 1/r
        w_: the current value(s) of du/dphi

    Returns:
        (du/dphi, dw/dphi)
    """
    return w_, -u_ + 1.5 * RS * (u_**2)


def temperature_to_rgb(t):
    """
    Convert a normalized disk temerature to an RGB color 
    (red -> orange -> yellow -> white)

    Parameters:
        t: array-like, temerature normalized to [0, 1], where 1.0 is the 
        hottest (inner disk) and 0.0 is the coolest (outerdisk).

    Returns:
        Array of shape (..., 3) with RGB values in [0, 1].
    """
    t = np.clip(t, 0.0, 1.0)
    stops = np.array([
        #[ t ,  R  ,  G  ,   B ]
        [0.00, 0.35, 0.02, 0.00],  # Dim deep red
        [0.50, 1.00, 0.60, 0.05],  # Orange
        [0.65, 1.00, 0.80, 0.20],  # Orange-yellow
        [0.85, 1.00, 0.90, 0.35],  # Yellow-white
        [0.98, 1.00, 1.00, 1.00],  # White
        [1.00, 0.80, 0.95, 1.00],  # Bright blue-white
    ])

    r = np.interp(t, stops[:, 0], stops[:, 1])
    g = np.interp(t, stops[:, 0], stops[:, 2])
    bl = np.interp(t, stops[:, 0], stops[:, 3])

    return np.stack([r, g, bl], axis=-1)


def apply_bloom(img, threshold, layers, boost):
    """
    Adds a glow around spesified pixels by blurring only the light above the 
    threshold and adding it ontop of the original image at several blur radii.

    Parameters:
        img: the image to be blurred.
        threshold: the threshold at which the pixels will be blurred.
        layers: the radii in which the pixels will be blurred around.
        boost: the factor the light is multiplied by bafore blurring.

    Returns:
        The new immage that has been blurred.
    """
    lum = img[..., 0] * 0.2126 + img[..., 1] * 0.7152 + img[..., 2] * 0.0722
    excess = np.clip(lum - threshold, 0.0, None) * boost
    bright_pass = img * (excess / np.maximum(lum, 1e-6))[..., None]

    glow = np.zeros_like(img)
    for sigma, weight in layers:
        glow += weight * gaussian_filter(bright_pass, sigma=(sigma, sigma, 0))

    return img + glow


# Setting up camera using provided constants
incl = np.radians(INCLINATION_DEG)
tilt = np.radians(TILT_DEG)
cam_pos = np.array([np.sin(incl), 0.0, np.cos(incl)]) * CAMERA_DIST

forward = -cam_pos / np.linalg.norm(cam_pos)
world_up = np.array([0.0, 0.0, 1.0])
if abs(np.dot(forward, world_up)) > 0.999:
    world_up = np.array([0.0, 1.0, 0.0])

right = np.cross(forward, world_up)
right /= np.linalg.norm(right)
up = np.cross(right, forward)

right, up = (np.cos(tilt) * right + np.sin(tilt) * up, 
             -np.sin(tilt) * right + np.cos(tilt) * up,)

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
b = r0 * sinpsi / np.sqrt(metric0)  # Impact parameter
w0 = u0 * np.sqrt(metric0) * (cospsi / sinpsi)  # du/dphi at phi=0

u = u0.copy()
w = w0.copy()
phi = np.zeros(N)

prev_z = pos0[:, 2].copy()
prev_pos = pos0.copy()

# Conditions to track for photon rays
active = np.ones(N, dtype=bool)
hit_disk = np.zeros(N, dtype=bool)
was_captured = np.zeros(N, dtype=bool)
result_disk_color = np.zeros((N, 3))
result_star_color = np.zeros((N, 3))
disk_hit_pos = np.zeros((N, 3))

# RK4 integration
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
    pos_new = r_new[:, None] * (cosphi[:, None] * e_r[idx] + sinphi[:, None]
                                * e_phi[idx])
    z_new = pos_new[:, 2]

    # Checks if the photon ray falls into the black hole or escapes into space
    captured = r_new <= RS * 1.001
    escaped = (r_new > ESCAPE_R)

    # Check if the photon ray impacts the accretion disk
    z_prev_sub = prev_z[idxs]
    crossed = ((np.sign(z_new) != np.sign(z_prev_sub)) & (r_new >= DISK_INNER)
                & (r_new <= DISK_OUTER))

    if crossed.any():
        t = z_prev_sub[crossed] / (z_prev_sub[crossed] - z_new[crossed])
        p_prev = prev_pos[idxs[crossed]]
        p_new = pos_new[crossed]

        hit_p = p_prev + t[:, None] * (p_new - p_prev)
        ci = idxs[crossed]
        hit_disk[ci] = True
        active[ci] = False
        disk_hit_pos[ci] = hit_p

    still = ~crossed
    if (captured & still).any():
        cidx = idxs[captured & still]
        active[cidx] = False
        was_captured[cidx] = True
    if (escaped & still & ~captured).any():
        eidx = idxs[escaped & still & ~captured]
        active[eidx] = False

    prev_z[idxs] = z_new
    prev_pos[idxs] = pos_new

# Color assignments
#result_disk_color[:] = 0.0
#result_star_color[:] = 0.0

# Creating background stars
bg_rand = rng.random(N)
is_star = bg_rand < STAR_DENSITY
star_brightness = rng.uniform(0.2, 1.0, size=N)
show_star = is_star & ~hit_disk & ~was_captured
result_star_color[show_star] = star_brightness[show_star, None]

# Disk shading
if hit_disk.any():
    hp = disk_hit_pos[hit_disk]
    r_hit = np.linalg.norm(hp, axis=1)
    r_hit = np.clip(r_hit, DISK_INNER, DISK_OUTER)

    # Temperature profile ~ r^(-3/4) (Shakura-Sunyaev-like folloff), 
    # normalized to [0, 1]
    temp_raw = r_hit ** (-0.75)
    tmin, tmax = (DISK_OUTER ** (-0.75)), (DISK_INNER ** (-0.75))
    temp_norm = (temp_raw - tmin) / (tmax - tmin)

    # Approximation of Doppler beaming
    v = np.sqrt(np.clip(RS / (2.0 * r_hit), 0.0, 0.98))
    gamma = 1.0 / np.sqrt(1.0 - (v**2))

    radial_hat = hp / r_hit[:, None]
    tang_hat = np.stack([-radial_hat[:, 1], radial_hat[:, 0], np.zeros_like(radial_hat[:, 0])], axis=-1)
    vel_hat = tang_hat  # Disk is assumed to rotate counter-clockwise viewed from +z

    los = cam_pos[None, :] - hp
    los = los / np.linalg.norm(los, axis=1, keepdims=True)

    cos_angle = np.einsum('ij,ij->i', vel_hat, los)
    doppler = 1.0 / (gamma * (1.0 - v*cos_angle))
    grav_redshift = np.sqrt(np.clip(1.0 - RS / r_hit, 1e-6, None))
    shift = doppler * grav_redshift

    brightness = np.clip(shift**3, 0.05, 8.0)
    brightness = brightness / (1.0 + 0.15 * brightness)

    color_shift = np.clip(0.5 + 5.0 * (shift - 1.0), 0.0, 1.0)
    rgb = temperature_to_rgb(np.clip(temp_norm + 0.25 * (color_shift - 0.5), 0, 1))

    disk_rgb = rgb * brightness[:, None]
    result_disk_color[hit_disk] = disk_rgb

# Creating the final image
#img = result_color.reshape(HEIGHT,WIDTH, 3)
star_img = result_star_color.reshape(HEIGHT,WIDTH, 3)
disk_img = result_disk_color.reshape(HEIGHT,WIDTH, 3)

if BLOOM_ENABLED:
    star_img = apply_bloom(star_img, STAR_BLOOM_THRESHOLD, STAR_BLOOM_LAYERS, STAR_BLOOM_BOOST)
    disk_img = apply_bloom(disk_img, DISK_BLOOM_THRESHOLD, DISK_BLOOM_LAYERS, DISK_BLOOM_BOOST)

img = star_img + disk_img
img = np.clip(img, 0.0, 1.0)
img = img ** (1.0/2.0)
img8 = (img * 255).astype(np.uint8)

out = Image.fromarray(img8, mode='RGB')
out.save('C:\\Users\\cedri\\Coding_Projects' + 
         '\\Black_Hole\\blackhole_lensing.png')
print('Saved C:\\Users\\cedri\\Coding_Projects' + 
      '\\Black_Hole\\blackhole_lensing.png')