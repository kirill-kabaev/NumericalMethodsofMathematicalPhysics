#!/usr/bin/env python
# coding: utf-8

# In[0]:


import matplotlib.pyplot as plt
import pylab
from matplotlib import mlab
import matplotlib as mpl
import numpy as np
#import math
from cmath import sqrt
from mpmath import *
get_ipython().run_line_magic('matplotlib', 'inline')

N = 1000
L=1
delta=1/N
V0 = 50
z=0.1

def V(x):
    return V0/((np.cosh((x-0.5)/z)) ** 2)

xlist = np.linspace(-10*L, 10*L, N)
Ulist = [V(x) for x in xlist]
xlist = np.asarray(xlist, dtype=float)
Ulist = np.asarray(Ulist, dtype=float)
d=0
def Tunel(E):
    global N
    global delta
    global d
    R=[0]*N
    f=[0]*N
    d=d+1
    #print(d)
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
plt.ylim(-1, V0+10)
plt.xlabel(u'x')
plt.ylabel(u'V(x)')


# In[1]:


import pylab 

V0=50
tun=[]
a=[]
b=[]
Emax=500
Emin=0.1
deltaE=(Emax-Emin)/N
Elist = np.linspace(Emin, Emax, N)
for i in range(N):
    tun.append(Tunel(Emin+deltaE*i))
    a.append(tun[i][0])
    b.append(tun[i][1])
pylab.plot (Elist, a, label = u'T(V0=50)')
#pylab.plot (Elist, b, label = u'R') 

V0 = 100
tun=[]
a=[]
b=[]
Emax=500
Emin=0.1
deltaE=(Emax-Emin)/N
Elist = np.linspace(Emin, Emax, N)
for i in range(N):
    tun.append(Tunel(Emin+deltaE*i))
    a.append(tun[i][0])
    b.append(tun[i][1])
pylab.plot (Elist, a, label = u'T(V0=100)')
#lab.plot (Elist, b, label = u'R') 

V0 = 150
tun=[]
a=[]
b=[]
Emax=500
Emin=0.1
deltaE=(Emax-Emin)/N
Elist = np.linspace(Emin, Emax, N)
for i in range(N):
    tun.append(Tunel(Emin+deltaE*i))
    a.append(tun[i][0])
    b.append(tun[i][1])
pylab.plot (Elist, a, label = u'T(V0=150)')
#lab.plot (Elist, b, label = u'R') 
 

pylab.legend(loc='lower right')
plt.xlabel(u'E')
plt.xlim(-1, Emax)


# In[2]:


def tLau(E):
    k=np.sqrt(E)
    return (np.sinh(np.pi*k*1/z)/np.sinh(np.pi*k*1/z)
pylab.plot (Elist,tlistLau, label = u'TLau')

