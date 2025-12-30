import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

# constants
c = 3e8
G = 6.67e-11
mass = 2e30 # mass can be adjusted and everything else will adjust itself

mass_ratio = (mass/(2e30))
rs = 2 * G * mass / (c**2)
dlambda = (rs/c) * 0.1

start = 9 # the number  of rs the photons will start from the blackhole
hight = 5 * rs
width = start/5 * hight

ph_num = 28
ph_range = np.linspace(hight, 0, ph_num)
# for a single orbiting photon use ph_range=np.linspace(2.5873539027501*rs, 2.6*rs, ph_num) 
# and mass=2e30 and start=9 and frames=300

class Blackhole:
    def __init__(self, x, y):
        self.x = x
        self.y = y

        ax.add_patch(plt.Circle((self.x, self.y), rs, color='r', zorder=2)) # event horizon
        ax.add_patch(plt.Circle((self.x, self.y), rs*1.5, color='c', fill = False)) # photon sphere
        ax.add_patch(plt.Circle((self.x, self.y), rs*2.6, color='c', fill = False)) # image
        

class Photon:
    def __init__(self, x, y):
        # cartesian coords
        self.pos = np.array([x, y])
        self.v   = np.array([-c, 0])
        
        # polar coords
        self.r       = np.hypot(x, y)
        self.theta   = np.atan2(y, x)
        self.dr      = (self.pos[0] * self.v[0] + self.pos[1] * self.v[1]) / self.r
        self.dtheta  = (self.pos[0] * self.v[1] - self.pos[1] * self.v[0]) / (self.r**2)
        self.ddr     = 0
        self.ddtheta = 0
        
        # trail history start
        self.path_x = [x]
        self.path_y = [y]

        self.trail, = ax.plot([], [], color='y', alpha=0.95, linewidth=1)
        self.photon = ax.add_patch(plt.Circle((self.pos[0], self.pos[1]), 100*mass_ratio, color='y'))

        # true if photon outside of bh
        self.active = True


    def ph_update(self):
        if (self.r < rs): 
            self.active = False
            return

        # convert update to cartesian
        self.pos[0] = np.cos(self.theta) * self.r
        self.pos[1] = np.sin(self.theta) * self.r

        # update trail
        self.path_x.append(self.pos[0])
        self.path_y.append(self.pos[1])
        
        self.photon.set_center((self.pos[0], self.pos[1]))
        self.trail.set_data(self.path_x, self.path_y)


    def geodesic(self, state):
        r, theta, dr, dtheta = state

        ddr = (r - 1.5 * rs) * (dtheta**2)
        ddtheta = -2 * dr * dtheta / r

        return np.array([dr, dtheta, ddr, ddtheta])


    def rk4(self):
        state = np.array([self.r, self.theta, self.dr, self.dtheta])

        k1 = self.geodesic(state)
        k2 = self.geodesic(state + dlambda/2*k1)
        k3 = self.geodesic(state + dlambda/2*k2)
        k4 = self.geodesic(state + dlambda*k3)

        self.r      += (dlambda/6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])
        self.theta  += (dlambda/6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
        self.dr     += (dlambda/6.0) * (k1[2] + 2*k2[2] + 2*k3[2] + k4[2])
        self.dtheta += (dlambda/6.0) * (k1[3] + 2*k2[3] + 2*k3[3] + k4[3])


# window setup
fig = plt.figure()
fig.set_facecolor('k')
fig.tight_layout()

ax = fig.add_subplot(aspect = 'equal')
ax.set_facecolor('k')
ax.set_xlim(-width/(start/5), width)
ax.set_ylim(-hight, hight)

ax.yaxis.set_label_position("right")
ax.yaxis.tick_right()
ax.tick_params(axis='both', colors='w') 

for spine in ax.spines.values():
    spine.set_color('w')

tick_xpos = np.arange(-width/(start/5), width+mass_ratio, rs)
tick_ypos = np.arange(-hight, hight+mass_ratio, rs)
ax.set_xticks(tick_xpos)
ax.set_yticks(tick_ypos)

tick_xlables = np.abs(np.round(tick_xpos/rs))
tick_ylables = np.abs(np.round(tick_ypos/rs))
ax.set_xticklabels(tick_xlables.astype(int))
ax.set_yticklabels(tick_ylables.astype(int))
ax.set_xlabel('Measured in Schwarzschild Radii', color='w')

# creating the objects
phs = [Photon(width-100*mass_ratio, i) for i in ph_range]
bh = Blackhole(0,0)
print('Finished creating photons and blackhole\nCreating animation...')

# animation
def animate(i):
    for ph in phs:
        if ph.active:
            ph.rk4()
            ph.ph_update()

    return [ph.trail for ph in phs]

ani = FuncAnimation(fig, animate, frames=300) 
#plt.show()
ani.save("blackhole.gif", writer=PillowWriter(fps=50), dpi=200)

print('Done')