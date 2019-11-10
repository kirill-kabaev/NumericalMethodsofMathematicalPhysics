#!/usr/bin/env python
# coding: utf-8

# In0]:
import matplotlib.pyplot as plt
import pylab
from matplotlib import mlab
import matplotlib as mpl
import numpy as np
from cmath import sqrt
from mpmath import *
get_ipython().run_line_magic('matplotlib', 'inline')

N=300
delta=1/N
dt=0.0001
k0=50
sigma=0.08
xc=0.3
fin=[]
finabsvkv0=[]
L=1
d=0.1*L
V0 = (1/(np.pi**(1/4)*(sigma)**(1/2)))**2
def V(x):
    if (abs((x-0.4*L))<=0.1*L/2) or (abs((x-0.4*L-2*0.1*L/2-d))<=0.1*L/2) : return V0
    else: return 0
for i in range(N):
    fin.append((1/(np.pi**(1/4)*(sigma)**(1/2)))*np.exp(-((delta*i-xc)**2)/(2*sigma**2) + 1j*k0*delta*i))                                        
def u(n):
    return -2-2*V(delta*n)*(delta**2) + 4*1j*(delta**2)/dt
               
xlist = np.linspace(-L, L, N)
UlistV = [V(x) for x in xlist]


for i in range(N):
    finabsvkv0.append((abs((1/(np.pi**(1/4)*(sigma)**(1/2)))*np.exp(-((delta*i-xc)**2)/(2*sigma**2) + 1j*k0*delta*i)))**2)

pylab.plot (xlist, UlistV)
pylab.plot (xlist, finabsvkv0)
plt.xlim(-L, L)
plt.ylim(0,10)
plt.xlabel(u'x')
plt.ylabel(u'V(x)')


# In[1]:
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
plt.style.use('seaborn-pastel')
fig = plt.figure()
ax = plt.axes(xlim=(-L, L), ylim=(0, 15))
line, = ax.plot([], [], lw=3)
 
def init():
    line.set_data([], [])
    return line,
def animate(j):
    x = np.linspace(-L, L, N)
    R=[0]*N
    S=[0]*N
    F=[0]*N
    fin[0]=0
    fin[N-1]=0
    for i in range(N-1):
        F[i]=-(fin[i+1]+fin[i-1]+np.conj(u(i))*fin[i])
    i=N
    while(i>=1):
        i=i-1
        R[i-1]=-1/(u(i)+R[i])
    i=N
    while(i>=1):
        i=i-1
        S[i-1]=-R[i-1]*(F[i]-S[i])
    fin[1]=S[0]
    for i in range(N-1):
        fin[i+1]=R[i]*fin[i]+S[i]
    evolution=[]
    for i in range(N):
        evolution.append((abs(fin[i]))**2)
    line.set_data(x, evolution)
    return line,
anim = FuncAnimation(fig, animate, init_func=init,
                               frames=250, interval=1, blit=True)
anim.save('wave2barier.gif', writer='imagemagick')
