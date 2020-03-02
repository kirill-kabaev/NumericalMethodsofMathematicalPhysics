#!/usr/bin/env python
# coding: utf-8

# In[1]:
import matplotlib.pyplot as plt
import pylab
from matplotlib import mlab
import matplotlib as mpl
import numpy as np
#import math
from cmath import sqrt
from mpmath import *
get_ipython().run_line_magic('matplotlib', 'inline')

p=0.1
N = 100
x_min = -6
x_max = 6
V0 = 20
Emax =V0
gamma = 0.5
a=1

def V(x):
    return -V0/((cosh(x/a)) ** 2)

xlist = np.linspace(x_min, x_max, N)
Ulist = [V(x) for x in xlist]

pylab.plot (xlist, Ulist)
pylab.plot (xlist, [-10]*N)
pylab.plot (xlist, [0]*N)
pylab.plot (xlist, [-20]*N)

e = 0.001
ex = 0.001

#функция вычисляет точки поворота
def PointReturn(E):
    x0=0
    points=[0,0]
    global x_min
    global x_max
    xmin = x_min
    xmax = x_max
    while((xmax-x0) >= ex):
        a = (xmax+x0)/2
        f1 = V(xmax)-E
        f2 = V(a)-E
        if(f1*f2 > 0):
            xmax = a
        else:
            x0 = a
    points[1] = x0
    
    while((x0-xmin) >= ex):
        a = (xmin+x0)/2
        f1 = V(xmin)-E
        f2 = V(a)-E
        if(f1*f2 > 0):
            xmin = a
        else:
            x0 = a
    points[0] = x0
    return(points)

def Integral(E):
    xpoints=[0,0]
    xpoints=PointReturn(E)
    S=0
    x=xpoints[0]
    nn=150
    h0=(xpoints[1]-xpoints[0])/nn
    def f(x_):
        return(sqrt(2*(E-V(x_))))
    while(x<xpoints[1]):
        S=S+(h0/3)*(f(x)+4*f(x+h0)+f(x+2*h0))
        x=x+2*h0
    S=S.real
    return(S)

n=0
c = 0.1
Emax = 19.9
En=[]
MainCondition = True
PRRes=[]
Sa0 = Integral(c)
Sb0 = Integral(Emax)
while(MainCondition):
    a0 = c
    b0 = Emax
    Funca0 = (Sa0/np.pi) - n - gamma
    Funcb0 = (Sb0/np.pi) - n - gamma
    if(Funca0*Funcb0>0):
        if(Funca0<0):
            break
    else:
        while((b0-a0)>=e):
            c=(a0+b0)/2
            Funca0 = (Integral(a0)/np.pi - n - gamma)
            Funcc = (Integral(c)/np.pi - n - gamma)
            if(Funca0*Funcc>0):
                a0=c
            else:
                b0=c     
    En.append(c)
    PRRes.append(PointReturn(c))
    n=n+1

plt.xlabel(u'x')
plt.ylabel(u'V(x)')

# In[2]:
import pylab 
print("n - En")
for i in range(len(En)):
    print(i,'-', En[i])
    xlist1 = np.linspace(PRRes[i][0], PRRes[i][1], N)
    pylab.plot (xlist1, [En[i]]*N, label = u'n='+str(i))
    pylab.legend(loc='lower right')
pylab.plot (xlist, Ulist)
plt.xlabel(u'x')
plt.ylabel(u'V(x)')

# In[3]:
figure = pylab.figure()
axes = figure.add_subplot (1, 1, 1)
Ulist0=[]
k=0
for k in range(n):
    Ulist0=[]
    xlist= np.linspace(PRRes[k][0], PRRes[k][1], N)
    xlist0 = np.concatenate([xlist, xlist]) 
    Ulist1 = [ Ulist0.append((2*(En[k] - V(x)))**(1/2)) for x in xlist]
    Ulist2 = [Ulist0.append(-(2*(En[k] - V(x)))**(1/2)) for x in xlist]
    pylab.plot (xlist0, Ulist0, label = u'n='+str(k))
pylab.legend(loc='upper right')
ax = plt.gca()
plt.xlabel(u'x')
plt.ylabel(u'p(x)')
plt.xlim(-7, 7)
plt.ylim(-7, 7)


