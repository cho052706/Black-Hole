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

ph_num = 10
ph_range = np.linspace(hight/2, 0, ph_num)

pull = True

class Blackhole:
    def __init__(self, x, y):
        self.x = x
        self.y = y

        ax.add_patch(plt.Circle((self.x, self.y), rs, color='r',))
        ax.add_patch(plt.Circle((self.x, self.y), rs*1.5, color='r', fill = False))
        ax.add_patch(plt.Circle((self.x, self.y), rs*2.6, color='r', fill = False))

    def pull(self, ph):
        if not ph.active: return

        force_dir = ph.pos - np.array([self.x, self.y]) 
        r = np.linalg.norm(force_dir)
        epsilon = 1
        fg = G * mass / (r**2 + epsilon)
        
        acceleration = -fg * (force_dir/r)
        v_new = ph.v + acceleration * dt
        ph.v = c * (v_new/(np.linalg.norm(v_new)))
        

class Photon:
    def __init__(self, x, y):
        self.pos = np.array([x, y], dtype=float)
        self.v = np.array([-c, 0], dtype=float) #Replace 30 with c at some point

        self.path_x = [x]
        self.path_y = [y]

        self.trail, = ax.plot([], [], color='white', alpha=1, linewidth=1)
        self.photon = ax.add_patch(plt.Circle((self.pos[0], self.pos[1]), 100, color='w'))

        self.active = True
        
        #self.pos = np.array(x, y)
        #self.v = np.array(-c, int(0))

    def draw(self):
        self.photon = ax.add_patch(plt.Circle((self.pos[0], self.pos[1]), 10, color='w'))
        self.trail = ax.add_line(Line2D((self.pos[0], self.pos[1]),((self.pos[0], self.pos[1])), color='w', lw=0.1))

    def ph_update (self):
        #if not self.active: return

        self.pos += self.v * dt
        
        self.r = np.hypot(self.pos[0], self.pos[1])
        if (self.r < rs): 
            self.active = False
            return

        self.path_x.append(self.pos[0])
        self.path_y.append(self.pos[1])
        
        self.photon.set_center((self.pos[0], self.pos[1]))
        self.trail.set_data(self.path_x, self.path_y)


#fig, ax = plt.subplots(constrained_layout = True)
fig = plt.figure()
fig.set_facecolor('k')
fig.tight_layout()

ax = fig.add_subplot(aspect = 'equal')
ax.set_facecolor('k')
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlim(-width, 2*width)                   # Edit (ax) size here
ax.set_ylim(-hight, hight)

phs = [Photon(2*width-100, i) for i in ph_range]

bh = Blackhole(0,0)

def animate(i):
    for ph in phs:
        if ph.active:
            ph.ph_update()
            bh.pull(ph)

    return [ph.trail for ph in phs]

ani = FuncAnimation(fig, animate, frames=300) 
ani.save("blackhole.gif", writer=PillowWriter(fps=50))
#plt.show()
print('done')