#!/usr/bin/env python
# coding: utf-8

# In[0]:


import matplotlib.pyplot as plt
import pylab
from matplotlib import mlab
import matplotlib as mpl
import numpy as np
from cmath import sqrt
from mpmath import *
get_ipython().run_line_magic('matplotlib', 'inline')

N = 1000
V0 = 100
L=1
delta=1/N

def V(x):
    if abs((x-L/2))<=0.8*L/2: return V0
    else: return 0

xlist = np.linspace(-L-L/3, L+L/3, N)
Ulist = [V(x) for x in xlist]
xlist = np.asarray(xlist, dtype=float)
Ulist = np.asarray(Ulist, dtype=float)

def Tunel(E):
    global N
    global delta
    global d
    R=[0]*N
    f=[0]*N
    f = np.asarray(f, dtype=complex)
    R = np.asarray(R, dtype=complex)
    def u(n):
        return -2 + ((delta)**2)*(E-V(delta*n))
    k=complex(sqrt(E))
    R[N-1]=-1/(u(N)/2+1j*k*delta)
    for i in reversed(range(N-1)):
        R[i]=-1/(u(i+1)+R[i+1])
    f[0] = (2j*k*delta/(R[0]+(u(0)/2 +1j*k*delta)));
    for i in range(N-1):
        f[i+1]=R[i]*f[i]
    return [(abs(f[N-1]))**2, (abs(f[0]-1))**2]

pylab.plot (xlist, Ulist)
plt.xlim(-1-L/3, L+L/3)
plt.ylim(-1, 120)
plt.xlabel(u'x')
plt.ylabel(u'V(x)')


# In[1]:


import pylab 
from ipywidgets import widgets
from IPython.display import display
tun=[]
a=[]
b=[]
Emin=0
Emax=500
deltaE=Emax/N
Elist = np.linspace(Emin, Emax, N)
nlist  = np.linspace(0, N-1, N)

for i in range(N):
    tun.append(Tunel(deltaE*i))
    a.append(tun[i][0])
    b.append(tun[i][1])
pylab.plot (Elist, a, label = u'T')   
plt.xlabel(u'E')
plt.xlim(-1, Emax)

def t(E):
    kk=np.sqrt(E)
    if E<V0: q=1j*np.sqrt(V0-E)
    else: q=np.sqrt(E-V0)
    k=kk/q
    Th=q*0.8*L
    return (abs(4*k/(((1+k)**2) * np.exp(-1j*Th)-((k-1)**2) * np.exp(1j*Th)))) ** 2
tlist = [t(E) for E in Elist]
pylab.plot (Elist,tlist, label = u'TSat')
pylab.legend(loc='lower right')
plt.axis('auto')



# In[2]:


diff=[]
for i in range(N):
    diff.append(abs(a[i]-tlist[i]))
print(np.amax(diff))
pylab.plot (nlist, diff)
plt.ylabel(u'diff')
plt.xlabel(u'iteration')
plt.xlim(100, 500)


# In[3]:


pylab.plot (Elist, a, label = u'T')
pylab.plot (Elist, b, label = u'R')    
plt.xlabel(u'E')
