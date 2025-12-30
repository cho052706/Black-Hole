import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.lines import Line2D

c = 3e8
G = 6.67e-11
dlambda = 9e-7

mass = 2e30
rs = 2 * G * mass / (c**2)

hight = 5 * rs
width = hight

ph_num = 20
ph_range = np.linspace(hight, 0, ph_num)

class Blackhole:
    def __init__(self, x, y):
        self.x = x
        self.y = y

        ax.add_patch(plt.Circle((self.x, self.y), rs, color='r', zorder=2)) # event horizon
        ax.add_patch(plt.Circle((self.x, self.y), rs*1.5, color='r', fill = False)) # photon sphere
        ax.add_patch(plt.Circle((self.x, self.y), rs*2.6, color='r', fill = False)) # image
        

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



    def ph_update(self):
        if (self.r < rs): 
            self.active = False
            return
        
        # update polar
        #self.r += self.dr * dt
        #self.theta += self.dtheta * dt

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

        #ddr = r * (dtheta**2) - ((c**2) * rs) / (2 * (r**2))
        #ddtheta = -2 * dr * dtheta / r

        #self.dr += self.ddr * dt
        #self.dtheta += self.ddtheta * dt

        return np.array([dr, dtheta, ddr, ddtheta])


    def rk4(self):
        state = np.array([self.r, self.theta, self.dr, self.dtheta])

        k1 = self.geodesic(state)
        k2 = self.geodesic(state + dlambda/2*k1)
        k3 = self.geodesic(state + dlambda/2*k2)
        k4 = self.geodesic(state + dlambda*k3)

        self.r += (dlambda/6.0) * (k1[0] + 2*k2[0] + 2*k3[0] + k4[0])
        self.theta += (dlambda/6.0) * (k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
        self.dr += (dlambda/6.0) * (k1[2] + 2*k2[2] + 2*k3[2] + k4[2])
        self.dtheta += (dlambda/6.0) * (k1[3] + 2*k2[3] + 2*k3[3] + k4[3])

        x = self.r*np.cos(self.theta)
        y = self.r*np.sin(self.theta)

        v = np.array([x, y])

        if (np.linalg.norm(v) > c):
            new_v = c * (v/np.linalg.norm(v))

            self.r = np.linalg.norm(new_v)
            self.theta = np.atan2(new_v[1], new_v[0])


fig = plt.figure()
fig.set_facecolor('k')
fig.tight_layout()

ax = fig.add_subplot(aspect = 'equal')
ax.set_facecolor('k')

ax.set_xlim(-width, width)
ax.set_ylim(-hight, hight)






phs = [Photon(width-100, i) for i in ph_range]
bh = Blackhole(0,0)

print('finished creating photons and blackhole')

def animate(i):
    for ph in phs:
        if ph.active:
            ph.rk4()
            ph.ph_update()

    return [ph.trail for ph in phs]


ani = FuncAnimation(fig, animate, frames=300) 
plt.show()

#ani.save("blackhole.gif", writer=PillowWriter(fps=50))

print('done')