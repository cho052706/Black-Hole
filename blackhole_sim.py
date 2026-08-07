import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter

# Constants
c = 3e8
G = 6.67e-11
mass = 2e30

mass_ratio = mass / (2e30)
rs = (2 * G * mass) / (c**2)
dlambda = 0.1 * (rs/c)

start = 9  # How far away the photon(s) start
height = 5 * rs
width = (start/5) * height

ph_num = 28
ph_range = np.linspace(height, 0, ph_num)
# For a single orbiting photon use:
# ph_range = [(((3*np.sqrt(3)) / 2) - (1.06e-2)) * rs] 


class Blackhole:
    """Draws a balck hole on the screen when initiated."""

    def __init__(self, x, y, ax, rs):
        """
        Generates one instance of a black hole using the provided (x, y)
        coordinates, plot, and Schwarzchild radius. A visual for the event 
        horizon, photon sphere, and the balck hole's shadow/image are also
        drawn.
        """
        self.x = x
        self.y = y
        self.rs = rs
 
        shadow_radius = ((3*np.sqrt(3)) / 2) * self.rs

        # Drawing the event horizon.
        ax.add_patch(plt.Circle((self.x, self.y), self.rs, color='r', 
                                zorder=2))
        # Drawing the photon sphere.
        ax.add_patch(plt.Circle((self.x, self.y), self.rs * 1.5, color='c', 
                                fill=False))
        # Drawing the shadow/image.
        ax.add_patch(plt.Circle((self.x, self.y), shadow_radius, color='c', 
                                fill=False, linestyle="--"))       


class Photon:
    """
    A simulated photon. This class contains methods that find the future 
    trajectory of the photon, inegrates the position of the photon, and 
    updates the position of the photon.
    """

    def __init__(self, x, y, ax, rs, dlambda, mass_ratio):
        """
        Generates one instance of a photon using the provided (x, y) position,
        plot, Schwarzschild radius, dlambda (used in rk4 integration), and the
        mass_ratio.
        """
        self.rs = rs
        self.dlambda = dlambda

        # Setting up cartesian coordinates.
        self.pos = np.array([x, y])
        self.v = np.array([-c, 0])
        
        # Setting up polar coordinates.
        self.r = np.hypot(x, y)
        self.theta = np.atan2(y, x)
        self.dr = ((self.pos[0]*self.v[0]) + (self.pos[1]*self.v[1])) / self.r
        self.dtheta = (((self.pos[0]*self.v[1]) - (self.pos[1]*self.v[0])) 
                       / (self.r**2))
        
        # Photon trail history start
        self.path_x = [x]
        self.path_y = [y]

        self.trail, = ax.plot([], [], color='y', alpha=0.95, linewidth=1)
        self.photon = ax.add_patch(plt.Circle((self.pos[0], self.pos[1]), 
                                              100 * mass_ratio, color='y'))

        # True if photon is outside of black hole, false otherwise.
        self.active = True

    def ph_update(self):
        """Updates the photon position"""
        if self.r < rs: 
            self.active = False
            return

        # Convert the update to cartesian coordinates.
        self.pos[0] = np.cos(self.theta) * self.r
        self.pos[1] = np.sin(self.theta) * self.r

        # Updates trail history.
        self.path_x.append(self.pos[0])
        self.path_y.append(self.pos[1])

        # Updates Matplotlib artists.
        self.photon.set_center((self.pos[0], self.pos[1]))
        self.trail.set_data(self.path_x, self.path_y)

    def geodesic(self, state):
        """
        Calcualtes future position based on Schwarzschild null geodesics
        and the current position/state of the photon.
        """
        r, theta, dr, dtheta = state

        ddr = (r - (1.5*rs)) * (dtheta**2)
        ddtheta = (-2 * dr * dtheta) / r

        return np.array([dr, dtheta, ddr, ddtheta])

    def rk4(self):
        """
        Integrates the position of the photon using Runge-Kutta 4 
        Integration.
        """
        state = np.array([self.r, self.theta, self.dr, self.dtheta])

        k1 = self.geodesic(state)
        k2 = self.geodesic(state + ((dlambda/2) * k1))
        k3 = self.geodesic(state + ((dlambda/2) * k2))
        k4 = self.geodesic(state + (dlambda*k3))

        self.r += (dlambda/6.0) * (k1[0] + (2*k2[0]) + (2*k3[0]) + k4[0])
        self.theta += (dlambda/6.0) * (k1[1] + (2*k2[1]) + (2*k3[1]) + k4[1])
        self.dr += (dlambda/6.0) * (k1[2] + (2*k2[2]) + (2*k3[2]) + k4[2])
        self.dtheta += (dlambda/6.0) * (k1[3] + (2*k2[3]) + (2*k3[3]) + k4[3])

# Window setup
fig = plt.figure()
fig.set_facecolor('k')
fig.tight_layout()

ax = fig.add_subplot(aspect = 'equal')
ax.set_facecolor('k')
ax.set_xlim(-width / (start/5), width)
ax.set_ylim(-height, height)

ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()
ax.tick_params(axis='both', colors='w') 

for spine in ax.spines.values():
    spine.set_color('w')

# Axis setup
tick_xpos = np.arange(-width / (start/5), width + mass_ratio, rs)
tick_ypos = np.arange(-height, height + mass_ratio, rs)
ax.set_xticks(tick_xpos)
ax.set_yticks(tick_ypos)

tick_xlables = np.abs(np.round(tick_xpos / rs))
tick_ylables = np.abs(np.round(tick_ypos / rs))
ax.set_xticklabels(tick_xlables.astype(int))
ax.set_yticklabels(tick_ylables.astype(int))
ax.set_xlabel('Measured in Schwarzschild Radii', color='w')

# Creating the balckhole and photons
phs = [Photon(width - (100*mass_ratio), i, ax, rs, dlambda, mass_ratio)
        for i in ph_range]
bh = Blackhole(0, 0, ax, rs)

print('Finished creating photons and blackhole\nCreating animation...')

# Creates the animation
def animate(i):
    artists = []
    for ph in phs:
        if ph.active:
            ph.rk4()
            ph.ph_update()

            artists.append(ph.photon)
            artists.append(ph.trail)

    return artists

ani = FuncAnimation(fig, animate, frames=300, blit=True, interval=20) 
ani.save("blackhole.gif", writer=PillowWriter(fps=50), dpi=200)

print('Done')