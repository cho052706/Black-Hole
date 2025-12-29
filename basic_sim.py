import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.lines import Line2D

c = 3*10**8
G = 6.67 * (10**-11)
dt = 10**-6

mass = 2*10**30
rs = 2 * G * mass / (c**2)

hight = 20000
width = 20000

ph_num = 20
ph_range = np.linspace(3/4*hight, -3/4*hight, ph_num)

class Blackhole:
    def __init__(self, x, y):
        self.x = x
        self.y = y

        ax.add_patch(plt.Circle((self.x, self.y), rs, color='r',))
        ax.add_patch(plt.Circle((self.x, self.y), rs*1.5, color='r', fill = False))
        ax.add_patch(plt.Circle((self.x, self.y), rs*2.6, color='r', fill = False))
        

class Photon:
    def __init__(self, x, y):
        # cartesian coords
        self.pos = np.array([x, y])
        self.v = np.array([-c, 0])
        
        # polar coords
        self.r = np.hypot(x, y)
        self.theta = np.atan2(y, x)
        self.dr = (self.pos[0] * self.v[0] + self.pos[1] * self.v[1]) / self.r
        self.dtheta = (self.pos[0] * self.v[1] - self.pos[1] * self.v[0]) / (self.r**2)
        self.ddr = 0
        self.ddtheta = 0
        
        # trail history start
        self.path_x = [x]
        self.path_y = [y]

        self.trail, = ax.plot([], [], color='white', alpha=1, linewidth=1)
        self.photon = ax.add_patch(plt.Circle((self.pos[0], self.pos[1]), 100, color='w'))

        # true if photon outside of bh
        self.active = True


    def ph_update (self):
        if (self.r < rs): 
            self.active = False
            return
        
        # update polar
        self.dr += self.ddr * dt
        self.dtheta += self.ddtheta * dt
        self.r += self.dr * dt
        self.theta += self.dtheta * dt

        # convert update to cartesian
        self.pos[0] = np.cos(self.theta) * self.r
        self.pos[1] = np.sin(self.theta) * self.r

        # update trail
        self.path_x.append(self.pos[0])
        self.path_y.append(self.pos[1])
        
        self.photon.set_center((self.pos[0], self.pos[1]))
        self.trail.set_data(self.path_x, self.path_y)


def geodesic (ph, rs=rs):
    r = ph.r
    theta = ph.theta
    dr = ph.dr
    dtheta = ph.dtheta

    ph.ddr = r * (dtheta**2) - ((c**2) * rs) / (2 * (r**2))
    ph.ddtheta = -2 * dr * dtheta / r


fig = plt.figure()
fig.set_facecolor('k')
fig.tight_layout()
ax = fig.add_subplot(aspect = 'equal')
ax.set_facecolor('k')
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlim(-width, 2*width)
ax.set_ylim(-hight, hight)

phs = [Photon(2*width-100, i) for i in ph_range]
bh = Blackhole(0,0)

print('finished creating photons and blackhole')

def animate(i):
    for ph in phs:
        if ph.active:
            geodesic(ph)
            ph.ph_update()

    return [ph.trail for ph in phs]


ani = FuncAnimation(fig, animate, frames=300) 
ani.save("blackhole.gif", writer=PillowWriter(fps=50))

print('done')