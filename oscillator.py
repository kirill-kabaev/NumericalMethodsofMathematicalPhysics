#In[1]:
import matplotlib.pyplot as plt
import pylab
from matplotlib import mlab
import matplotlib as mpl
import numpy as np
from numpy import (array,zeros,ones,arange,linspace,logspace,
                   float64,int64,sin,cos,pi,exp,log,sqrt,abs,
                   nan,inf,any,all,sort,hstack,vstack,hsplit,
                   delete,insert,append,eye,fromfunction,
                   trace,diag,average,std,outer,meshgrid)
from numpy.linalg import det,inv,solve,eig
from cmath import sqrt
from mpmath import *
get_ipython().run_line_magic('matplotlib', 'inline')
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
plt.style.use('seaborn-pastel')
from scipy.sparse import dia_matrix
from scipy.sparse import diags
from scipy import sparse
from scipy.sparse.linalg import inv
from scipy.sparse import csc_matrix

def startconditiongauss():
    global jlist, Elist, Sigmalist1, Sigmalist, Sigmalist2, Psi0absvkv,Psi0, I, E, Sigma, N, dt, Psi0absvkvnorm, tlist, T
    N=100
    Psi0=[]
    dt=0.01
    a=round(N/10)

    jlist = np.linspace(1, N, N)
    jlist = np.asarray(jlist, dtype=float)
    Psi0=np.exp(-((jlist-round(N/2))**2)/(2*a**2))
    Psi0absvkv=abs(Psi0*np.conj(Psi0))
    S0=0
    j1=1
    jN=len(jlist)
    h0=1
    T=10000

    tlist=np.linspace(1, T, T)
    tlist = np.asarray(tlist, dtype=float)

    while(j1<jN-1): 
        S0=S0+(h0/3)*(Psi0absvkv[j1]+4*Psi0absvkv[j1+h0]+Psi0absvkv[j1+2*h0])
        j1=j1+2*h0
    Psi0norm=np.sqrt((1/S0))*Psi0
    Psi0absvkvnorm=(1/S0)*Psi0absvkv


startconditiongauss()
fig, (ax1) = plt.subplots(
    nrows=1, ncols=1,
    figsize=(10, 4)
)
ax1.scatter(x=jlist, y=Psi0absvkvnorm, marker='o', c='b')
ax1.set_xlabel('$j$')
ax1.set_ylabel('$|Psi(t=0)|^2$')
ax1.set_ylim(0, max(Psi0absvkvnorm)+ max(Psi0absvkvnorm)/10)

#In[2]:
#MathWaiting
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
plt.style.use('seaborn-pastel')
from scipy.sparse import dia_matrix
from scipy.sparse import diags
from scipy import sparse
from scipy.sparse.linalg import inv
from scipy.sparse import csc_matrix
startconditiongauss()
p=0

fig = plt.figure()
ax = plt.axes(xlim=(0, len(tlist)), ylim=(0,  100))
scatter, = ax.plot([], [],'bo', lw=3)

def init():
    scatter.set_data([], [])
    return scatter,
def animate(j):
    global p, Psi0, Psi0absvkv, tlist
    p=p+1
    s=0
    Ww=10
    p=p+1 
    s=1
    w=1
    A=10
    Hminus1=array([0.]*(N-1))
    Hplus1=array([0.]*(N-1))
    H0=array([0.]*(N))
    for i in range(N):
        s=s+1
        H0[i]=Ww*s
        if i!=N-1:
            Hminus1[i]=np.sqrt(s)*A*np.cos((w)*(p-1)*dt)
            Hplus1[i]=np.sqrt(s)*A*np.cos((w)*(p-1)*dt)
    I = csc_matrix(sparse.eye(N).toarray())
    diagonals = [H0, Hminus1, Hplus1]
    H=csc_matrix(diags(diagonals, [0, -1, 1]).toarray())
    U=(I-(1j*H*dt)/2)*inv((I+(1j*H*dt)/2))
    S0=0
    j1=1
    jN=len(jlist)
    h0=1
    while(j1<jN-1): 
        S0=S0+(h0/3)*(Psi0absvkv[j1]+4*Psi0absvkv[j1+h0]+Psi0absvkv[j1+2*h0])
        j1=j1+2*h0
    S0=abs(S0)
    Psi0absvkvnorm=(1/S0)*abs(Psi0*np.conj(Psi0))        
    E=jlist@Psi0absvkvnorm 
    t=p-1
    Elist[t]=E
    scatter.set_data(tlist, Elist)
    Psit=U.dot(Psi0)
    Psi0=Psit 
    return scatter,
anim = FuncAnimation(fig, animate, init_func=init,
                               frames=T, interval=2, blit=True)
anim.save('oscillevolutemathwaitng.gif', writer='imagemagick')

#In[3]:
#Sigma
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
plt.style.use('seaborn-pastel')
from scipy.sparse import dia_matrix
from scipy.sparse import diags
from scipy import sparse
from scipy.sparse.linalg import inv
from scipy.sparse import csc_matrix
startconditiongauss()
p=0


fig = plt.figure()
ax = plt.axes(xlim=(0, len(tlist)), ylim=(0, 60))
scatter, = ax.plot([], [],'bo', lw=3)
def init():
    scatter.set_data([], [])
    return scatter,
def animate(j):
    global p, Psi0, Psi0absvkv, Sigmalist
    p=p+1
    s=0
    Ww=10
    p=p+1 
    s=1
    w=1
    A=10
    Hminus1=array([0.]*(N-1))
    Hplus1=array([0.]*(N-1))
    H0=array([0.]*(N))

    for i in range(N):
        s=s+1
        H0[i]=Ww*s
        if i!=N-1:
            Hminus1[i]=np.sqrt(s)*A*np.cos((w)*(p-1)*dt)
            Hplus1[i]=np.sqrt(s)*A*np.cos((w)*(p-1)*dt)
    I = csc_matrix(sparse.eye(N).toarray())
    diagonals = [H0, Hminus1, Hplus1]
    H=csc_matrix(diags(diagonals, [0, -1, 1]).toarray())
    U=(I-(1j*H*dt)/2)*inv((I+(1j*H*dt)/2))
    S0=0
    j1=1
    jN=len(jlist)
    h0=1
    while(j1<jN-1): 
        S0=S0+(h0/3)*(Psi0absvkv[j1]+4*Psi0absvkv[j1+h0]+Psi0absvkv[j1+2*h0])
        j1=j1+2*h0
    S0=abs(S0)
    #print(S0)
    Psi0absvkvnorm=(1/S0)*abs(Psi0*np.conj(Psi0))        
    E=jlist@Psi0absvkvnorm
    Sigma=sum(((jlist-E)**2)/len(jlist))
    
    t=p-1
    Sigmalist[t]=np.sqrt(Sigma)
    scatter.set_data(tlist, Sigmalist)
    Psit=U.dot(Psi0)
    Psi0=Psit 
    return scatter,
anim = FuncAnimation(fig, animate, init_func=init,
                               frames=T, interval=1, blit=True)
anim.save('oscillevoluteSigma.gif', writer='imagemagick')

#In[4]:
#\psi\^2
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
plt.style.use('seaborn-pastel')
from scipy.sparse import dia_matrix
from scipy.sparse import diags
from scipy import sparse
from scipy.sparse.linalg import inv
from scipy.sparse import csc_matrix


startconditiongauss()
p=0
fig = plt.figure()
ax = plt.axes(xlim=(0, len(jlist)), ylim=(0,  0.2))
scatter, = ax.plot([], [],'bo', lw=3)
def init():
    scatter.set_data([], [])
    return scatter,
def animate(j):
    global Psi0,Psi0absvkv
    global p
    p=p+1
    s=0
    Ww=1
    w=10
    A=10
    
    Hminus1=array([0.]*(N-1))
    Hplus1=array([0.]*(N-1))
    H0=array([0.]*(N))

    for i in range(N):
        s=s+1
        H0[i]=Ww*s
        if i!=N-1:
            Hminus1[i]=np.sqrt(s)*A*np.cos((w)*(p-1)*dt)
            Hplus1[i]=np.sqrt(s)*A*np.cos((w)*(p-1)*dt)
    I = csc_matrix(sparse.eye(N).toarray())
    diagonals = [H0, Hminus1, Hplus1]
    H=csc_matrix(diags(diagonals, [0, -1, 1]).toarray())
    U=(I-(1j*H*dt)/2)*inv((I+(1j*H*dt)/2))
    S0=0
    j1=1
    jN=len(jlist)
    h0=1
    while(j1<jN-1): 
        S0=S0+(h0/3)*(Psi0absvkv[j1]+4*Psi0absvkv[j1+h0]+Psi0absvkv[j1+2*h0])
        j1=j1+2*h0
    S0=abs(S0)
    #print(S0)
    Psi0absvkvnorm=(1/S0)*abs(Psi0*np.conj(Psi0))
    scatter.set_data(jlist, Psi0absvkvnorm)
    Psit=U.dot(Psi0)
    Psi0=Psit  
    Psi0absvkv=abs(Psi0*np.conj(Psi0))
    
    #scatter.set_data(jlist, Psi0norm)
    return scatter,
anim = FuncAnimation(fig, animate, init_func=init,
                               frames=T, interval=1, blit=True)
anim.save('oscillevolutepsi2.gif', writer='imagemagick')

